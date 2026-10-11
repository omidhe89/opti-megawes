% MODEL_FOR_OPTIMIZATION Defines the dynamic and algebraic constraints for a single time step in the OCP.
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 14, 2025
% Last Modified: June 14, 2025

function [Pk, omega_dot, v_dot, lift_force, winch_acc_scaled,  alpha, beta, g, lbg, ubg, equality] = model_for_optimization(symbolic_funcs, Xk, X0, Uk, Zk, Z0, Phi, g, lbg, ubg, params_num, paramsColumnVec, ocp_bounds, equality, h)

%
% This function sets up the constraints for the nonlinear programming problem.
% It includes model dynamics (ODE), algebraic constraints, and various path and operational limits.
% The ODE defects use the implicit MIDPOINT rule (2nd order, time-symmetric,
% no artificial numerical dissipation -- important because the objective
% integrates mechanical power, which plain implicit Euler biases at O(h)).
%
% Inputs:
%   symbolic_funcs  - Structure containing symbolic CasADi functions for model components.
%   Xk              - CasADi symbolic variable for states at time step k+1.
%   X0              - CasADi symbolic variable for states at time step k.
%   Uk              - CasADi symbolic variable for control inputs at time step k.
%   Zk              - CasADi symbolic variable for algebraic states at time step k+1.
%   Z0              - CasADi symbolic variable for algebraic states at time step k
%                     (needed to evaluate the tether force at the interval midpoint).
%   Phi             - CasADi symbolic variable for homotopy parameters.
%   g               - Cell array of CasADi symbolic expressions for constraints (updated).
%   lbg             - Numeric array of lower bounds for constraints (updated).
%   ubg             - Numeric array of upper bounds for constraints (updated).
%   params_num      - Structure with numerical parameters of the system.
%   paramsColumnVec - Column vector of numerical parameters for CasADi functions.
%   ocp_bounds      - Structure with bounds for OCP variables.
%   equality        - Boolean array indicating equality constraints (updated).
%   h               - Step size for the current interval.
%
% Outputs:
%   Pk              - Scaled power output at the right node (path constraint
%                     and energy/flatness quadrature in building_NLP.m).
%   omega_dot       - Scaled angular accelerations.
%   lift_force      - Calculated lift force.
%   winch_acc_scaled- Scaled winch acceleration.
%   alpha           - Angle of attack at the interval midpoint (SI units).
%   beta            - Sideslip angle at the interval midpoint (SI units).
%   g               - Updated cell array of CasADi symbolic expressions for constraints.
%   lbg             - Updated numeric array of lower bounds for constraints.
%   ubg             - Updated numeric array of upper bounds for constraints.
%   equality        - Updated boolean array indicating equality constraints.

    % Extract environmental parameters
    windDir = params_num.env.wind.direction;
    baseWindSpeed = params_num.env.wind.base_velocity;

    % Extract dimensions
    Nx = params_num.Nx;    % Number of states (x(1:2) -> winch, x(3:15) -> kite)
    Nz = params_num.Nz;    % Number of algebraic states (e.g., [Magnitude, el, az])
    Nu = params_num.Nu;    % Number of control inputs (according to the chosen control mode)

    % Extract OCP bounds for various quantities
    alpha_min = ocp_bounds.alpha_min;
    alpha_max = ocp_bounds.alpha_max;
    beta_min = ocp_bounds.beta_min;
    beta_max = ocp_bounds.beta_max;
    v_app_min = ocp_bounds.v_app_min;
    v_app_max = ocp_bounds.v_app_max;
    winch_acc_min = ocp_bounds.winch_acc_min;
    winch_acc_max = ocp_bounds.winch_acc_max;
    kite_acc_min = ocp_bounds.kite_acc_min;
    kite_acc_max = ocp_bounds.kite_acc_max;
    power_min = ocp_bounds.power_min;
    power_max = ocp_bounds.power_max;

    x_L = ocp_bounds.x_L;
    x_U = ocp_bounds.x_U;
    u_L = ocp_bounds.u_L;
    u_U = ocp_bounds.u_U;
    z_L = ocp_bounds.z_L;
    z_U = ocp_bounds.z_U;

    % Convert scaled decision variables to SI units
    X0_si = scaled_to_si(X0, x_L, x_U);
    Uk_si = scaled_to_si(Uk, u_L, u_U);
    Zk_si = scaled_to_si(Zk, z_L, z_U);
    Z0_si = scaled_to_si(Z0, z_L, z_U);
    Xk_si = scaled_to_si(Xk, x_L, x_U);

    % Interval-midpoint state and algebraic values for the implicit-midpoint
    % integrator. z at the midpoint is the average of the two node values --
    % using only Zk here would degrade the scheme back to 1st order through
    % the tether-force coupling.
    X_mid_si = 0.5 .* (X0_si + Xk_si);
    Z_mid_si = 0.5 .* (Z0_si + Zk_si);

    % Get kite data and tether dynamics
    % M_OB: rotation matrix from body frame to world frame
    % (pos_W is not needed here: the position states x(12:14) ARE the
    % wind-frame Cartesian position, and the flight-region limits are plain
    % box bounds on those decision variables -- see system_OCP_bounds.m. The
    % former pos_W path constraints would be exactly redundant with them.)
    [~, ~, M_OB] = symbolic_funcs.kite_data(Xk_si, paramsColumnVec);
    % Algebraic (tether) constraint is enforced ONCE per node: the right node
    % of each interval here, node 0 in building_NLP.m. Do NOT also enforce it
    % at the left node -- every interior node is the right node of one
    % interval and the left node of the next, so that duplicates the rows and
    % makes the constraint Jacobian rank-deficient (LICQ violation). With
    % alg_eq = 0 holding at both nodes, the averaged Z_mid below matches the
    % exact midpoint algebraic value to O(h^2), which is all the 2nd-order
    % midpoint scheme needs.
    [alg_eq, T_kite, ~] = symbolic_funcs.tether_dynamics(Xk_si, Zk_si, paramsColumnVec);
    tether_forces.kite = T_kite;
    tether_forces.winch = Zk_si(1);

    % Tether-fuselage contact avoidance constraints
    if params_num.sim.ablation_study == 0
        % Transform kite tether force from world frame to body frame
        T_kite_O = transformFromWtoO(windDir,T_kite);
        T_kite_B = M_OB'*T_kite_O; % Rotate from inertial (O) to body (B) frame

        % Calculate angles (theta_kite, phi_kite) to avoid tether contact
        theta_kite = atan2(T_kite_B(3), T_kite_B(1)); % Angle in x-z plane
        g = {g{:}, si_to_scaled(theta_kite, 0, pi)};
        lbg = [lbg; si_to_scaled(deg2rad(9.5), 0, pi)]; % 2.5 degrees
        ubg = [ubg; si_to_scaled(deg2rad(170.5), 0, pi)]; % 177.5 degrees
        equality = [equality; false(1, 1)];

        phi_kite = atan2(T_kite_B(3), T_kite_B(2)); % Angle in y-z plane
        g = {g{:}, si_to_scaled(phi_kite, 0, pi)};
        lbg = [lbg; si_to_scaled(deg2rad(15), 0, pi)]; % 2.5 degrees
        ubg = [ubg; si_to_scaled(deg2rad(165), 0, pi)]; % 177.5 degrees
        equality = [equality; false(1, 1)];
    end

    % Calculate aerodynamic variables (apparent wind speed, angle of attack,
    % sideslip angle) at the interval midpoint, where the aerodynamic forces
    % enter the dynamics.
    [Va, alpha, beta] = symbolic_funcs.aerodynamics_vars(X_mid_si, paramsColumnVec);

    % Constraints on aerodynamic variables
    g = {g{:}, si_to_scaled(Va, v_app_min, v_app_max)};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    g = {g{:}, si_to_scaled(alpha, alpha_min, alpha_max)};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];
    g = {g{:}, si_to_scaled(beta, beta_min, beta_max)};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Constraint on generated power (path bound at the right NODE).
    % NOTE: keep the bound, the energy quadrature and the flatness penalty
    % (building_NLP.m) all on NODE power. Midpoint power r*Z_mid(1)*X_mid(2)
    % is built from averages of adjacent nodes, so an odd-even (Nyquist) node
    % oscillation cancels out of it -- with any of these terms evaluated at
    % the midpoint the optimizer exploits that blind spot (bang-bang winch
    % acceleration, sawtooth node power invisible to the objective).
    Pk_si = params_num.winch.radius * Zk_si(1) * Xk_si(2); % Power in SI units
    Pk = si_to_scaled(Pk_si, power_min, power_max); % Scale power
    g = {g{:}, vertcat(Pk)};
    lbg = [lbg; 0]; % Lower bound (scaled)
    ubg = [ubg; 1]; % Upper bound (scaled)
    equality = [equality; false(1, 1)];

    % Control mode-specific dynamics and constraints
    if params_num.ctrl_mode == 1 % Aerodynamic force and torque control
        M_winch = Uk_si(1);      % Winch moment
        F_aero_B = Uk_si(2:4);   % Aerodynamic force in body frame (control variable)
        M_aero_B = Uk_si(5:7);   % Aerodynamic moment in body frame (control variable)
        cs_defl = Uk_si(8:10);   % Control surface deflections
        F_prop_B = Uk_si(end);   % Propeller force

        % Enforce that F_aero_B and M_aero_B come from the aerodynamic model.
        % Evaluated at the MIDPOINT state: these forces enter the dynamics rhs
        % at X_mid, so tying them to the aero model at the right node would
        % inject an O(h) stage error and degrade the midpoint scheme to 1st
        % order. (cs_defl is a per-interval control, i.e. its stage value.)
        [F_a_B_expr,  M_a_B_expr] = symbolic_funcs.F_M_aero(X_mid_si, cs_defl, paramsColumnVec);
        g = {g{:}, (F_aero_B - F_a_B_expr)./ abs((u_U(2:4) - u_L(2:4))/2)}; % Force residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];
        g = {g{:}, (M_aero_B - M_a_B_expr) ./ abs((u_U(5:7) - u_L(5:7))/2)}; % Moment residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        % Implicit-midpoint integration for ODE (dynamics at the averaged state)
        rhs = vertcat(symbolic_funcs.winch_dynamics(X_mid_si, M_winch, Z_mid_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(X_mid_si, Z_mid_si, F_aero_B, M_aero_B, F_prop_B, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    elseif params_num.ctrl_mode == 2 % Fictitious force and torque control
        M_winch = Uk_si(1);      % Winch moment
        F_fict_B = Uk_si(2:4);   % Fictitious force in body frame (control variable)
        M_fict_B = Uk_si(5:7);   % Fictitious moment in body frame (control variable)
        F_aero_B = zeros(3,1);   % Aerodynamic force is zero in this mode
        M_aero_B = zeros(3,1);   % Aerodynamic moment is zero in this mode

        % Implicit-midpoint integration for ODE with fictitious forces
        rhs = vertcat(symbolic_funcs.winch_dynamics(X_mid_si, M_winch, Z_mid_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(X_mid_si, Z_mid_si, F_fict_B, M_fict_B, 0, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    else % Homotopy control (blend of aerodynamic and fictitious forces)
        M_winch = Uk_si(1);      % Winch moment
        F_fic_a_B = Uk_si(2:4);  % Fictitious/aerodynamic force (depending on Phi)
        M_fic_a_B = Uk_si(5:7);  % Fictitious/aerodynamic moment (depending on Phi)
        cs_defl = Uk_si(8:10);   % Control surface deflections
        F_prop_B = Uk_si(11);    % Propeller force
        F_aero_B = Uk_si(12:14); % Actual aerodynamic force (used for consistency constraints)
        M_aero_B = Uk_si(15:17); % Actual aerodynamic moment (used for consistency constraints)

        % Enforce that F_aero_B and M_aero_B come from the aerodynamic model.
        % Evaluated at the MIDPOINT state -- same reasoning as ctrl_mode 1:
        % the blended force enters the rhs at X_mid.
        [F_a_B_expr,  M_a_B_expr] = symbolic_funcs.F_M_aero(X_mid_si, cs_defl, paramsColumnVec);
        g = {g{:}, (F_aero_B - F_a_B_expr) ./ abs((u_U(12:14) - u_L(12:14))/2)}; % Force residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        g = {g{:}, (M_aero_B - M_a_B_expr) ./ abs((u_U(15:17) - u_L(15:17))/2)}; % Moment residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        % Combine fictitious/aerodynamic forces/moments based on homotopy parameter Phi(1)
        F_M_hat_a_B = Phi(1) .* [F_fic_a_B; M_fic_a_B] + (1 - Phi(1)) .* [F_aero_B; M_aero_B];
        % Apply homotopy to propeller force based on Phi(2)
        F_hat_p_B = Phi(2) .* F_prop_B;

        % Implicit-midpoint integration for ODE with blended forces
        rhs = vertcat(symbolic_funcs.winch_dynamics(X_mid_si, M_winch, Z_mid_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(X_mid_si, Z_mid_si, F_M_hat_a_B(1:3), F_M_hat_a_B(4:6), F_hat_p_B, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    end

    % NOTE on the quaternion rows (x9:11 = [q1;q2;q3], x15 = q0): under the
    % implicit MIDPOINT rule they need no special treatment -- evaluating the
    % kinematics q_dot = 0.5*Om(w)*q at the averaged state makes those defect
    % rows read
    %   q_next - q_prev - (h/4)*Om(w_mid)*(q_prev + q_next) = 0,
    % which is the Cayley scheme: dotting with (q_prev + q_next) kills the
    % skew-symmetric Om term identically, so ||q_next|| = ||q_prev|| EXACTLY
    % for any body rate and any step size. A single unit-norm constraint at
    % node 0 (building_NLP.m) therefore pins ||q|| = 1 at every node. Under a
    % plain implicit-EULER defect this is impossible: the step is orthogonal
    % to q_next, so unit norm at two consecutive nodes forces w = 0 (no
    % rotation) and the NLP is infeasible for any flying trajectory.

    % Constraint on translational acceleration magnitude (midpoint accelerations)
    trans_acc = rhs(3:5)' * rhs(3:5); % Squared magnitude of translational acceleration
    g = {g{:}, si_to_scaled(trans_acc, kite_acc_min, (kite_acc_max * params_num.env.g)^2)}; % Scaled constraint
    % No lower bound: trans_acc is a square, so a bound at 0 is vacuous, but IPOPT
    % still puts a barrier on it, which degenerates wherever rhs(3:5) ~ 0 (e.g.
    % steady circling, where the body-frame velocity is nearly constant).
    lbg = [lbg; -inf];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Constraint on winch acceleration
    winch_acc = rhs(2); % Winch acceleration from RHS
    winch_acc_scaled = si_to_scaled(winch_acc , winch_acc_min, winch_acc_max); % Scaled winch acceleration
    g = {g{:}, winch_acc_scaled};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Combine algebraic and ordinary differential equation (ODE) constraints
    % These are equality constraints enforcing model dynamics.
    % NOTE on the alg_eq scaling: the residual is eq = p - p0, i.e. all three
    % components are POSITION mismatches in meters (W-frame x,y,z; see
    % tether_quasi_static_single_shooting.m). They correspond axis-for-axis to
    % the Cartesian position states x(12:14) (same W-frame axes, same units),
    % so scale by those state ranges -- the same convention as the neighboring
    % ode_eq ./ (x_U - x_L) term (and the node-0 copy in building_NLP.m).
    alg_scale = x_U(12:14) - x_L(12:14); % per-axis flight-region extent (m)
    g = {g{:}, vertcat(alg_eq ./ alg_scale, ode_eq ./ (x_U - x_L))};
    lbg = [lbg; zeros(Nx+Nz,1)];
    ubg = [ubg; zeros(Nx+Nz,1)];
    equality = [equality; true(Nx+Nz, 1)];

    % NOTE: deliberately NO per-node unit-norm constraint on the quaternion.
    % The midpoint defect propagates ||q|| EXACTLY from node to node (see the
    % note above), so the single node-0 constraint suffices. Re-adding a norm
    % constraint per node would be exactly redundant with the quaternion
    % defect rows and make the constraint Jacobian rank-deficient (LICQ
    % violation -> exploding dual variables, unbounded Hessian regularization,
    % MUMPS failures).

    % Calculate and scale angular velocity derivatives and lift force
    v_dot = si_to_scaled(rhs(3:5), -4*params_num.env.g, 4*params_num.env.g); % range the v_dot cost weight is tuned for (as in awes_dms)
    omega_dot = si_to_scaled(rhs(6:8), -0.5 * ones(3,1), 0.5 * ones(3,1)); % Scaled angular accelerations

    lift_force = -F_aero_B(1)*sin(alpha) - F_aero_B(3)*cos(alpha); % Lift force calculation (midpoint stage)

end
