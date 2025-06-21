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

function [Pk, omega_dot, lift_force, winch_acc_scaled,  alpha, beta, g, lbg, ubg, equality] = model_for_optimization(symbolic_funcs, Xk, X0, Uk, Zk, Phi, g, lbg, ubg, params_num, paramsColumnVec, ocp_bounds, equality, h)

%
% This function sets up the constraints for the nonlinear programming problem.
% It includes model dynamics (ODE), algebraic constraints, and various path and operational limits.
%
% Inputs:
%   symbolic_funcs  - Structure containing symbolic CasADi functions for model components.
%   Xk              - CasADi symbolic variable for states at time step k+1.
%   X0              - CasADi symbolic variable for states at time step k.
%   Uk              - CasADi symbolic variable for control inputs at time step k.
%   Zk              - CasADi symbolic variable for algebraic states at time step k+1.
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
%   Pk              - Scaled power output.
%   omega_dot       - Scaled angular accelerations.
%   lift_force      - Calculated lift force.
%   winch_acc_scaled- Scaled winch acceleration.
%   alpha           - Angle of attack (SI units).
%   beta            - Sideslip angle (SI units).
%   g               - Updated cell array of CasADi symbolic expressions for constraints.
%   lbg             - Updated numeric array of lower bounds for constraints.
%   ubg             - Updated numeric array of upper bounds for constraints.
%   equality        - Updated boolean array indicating equality constraints.

    % Extract environmental parameters
    windDir = params_num.env.wind.direction;
    baseWindSpeed = params_num.env.wind.base_velocity;

    % Extract dimensions
    Nx = params_num.Nx;    % Number of states (e.g., x(1:2) -> winch, x(3:14) -> ap4)
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
    pos_W_x_min = ocp_bounds.pos_W_x_min;
    pos_W_x_max = ocp_bounds.pos_W_x_max;
    pos_W_y_min = ocp_bounds.pos_W_y_min;
    pos_W_y_max = ocp_bounds.pos_W_y_max;
    pos_W_z_min = ocp_bounds.pos_W_z_min;
    pos_W_z_max = ocp_bounds.pos_W_z_max;
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
    Xk_si = scaled_to_si(Xk, x_L, x_U);

    % Get kite data and tether dynamics
    % pos_W: position of kite in wind frame
    % M_OB: rotation matrix from body frame to world frame
    [~, pos_W, ~, M_OB, ~] = symbolic_funcs.kite_data(Xk_si, paramsColumnVec);
    [alg_eq, T_kite, ~] = symbolic_funcs.tether_dynamics(Xk_si, Zk_si, paramsColumnVec);
    tether_forces.kite = T_kite;
    tether_forces.winch = Zk_si(1);

    % Constraints on kite position in the world frame
   
    % Note: The first position constraint on pos_W(1) is commented out.
    % g = {g{:}, si_to_scaled(pos_W(1), pos_W_x_min, pos_W_x_max)};
    % lbg = [lbg; 0];
    % ubg = [ubg; 1];
    % equality = [equality; false(1, 1)]; 
                                     

    % Constraint on the y-position of the kite based on path type
    g = {g{:}, si_to_scaled(pos_W(2), pos_W_y_min , pos_W_y_max)};
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Constraint on the z-position (altitude) of the kite
    g = {g{:}, si_to_scaled(pos_W(3), pos_W_z_min, pos_W_z_max)}; %160
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Tether-fuselage contact avoidance constraints
    if params_num.sim.ablation_study == 0
        % Transform kite tether force from world frame to body frame
        T_kite_O = transformFromWtoO(windDir,T_kite);
        T_kite_B = M_OB'*T_kite_O; % Rotate from inertial (O) to body (B) frame

        % Calculate angles (theta_kite, phi_kite) to avoid tether contact
        theta_kite = atan2(T_kite_B(3), T_kite_B(1)); % Angle in x-z plane
        g = {g{:}, si_to_scaled(theta_kite, 0, pi)};
        lbg = [lbg; si_to_scaled(deg2rad(5), 0, pi)]; % 5 degrees
        ubg = [ubg; si_to_scaled(deg2rad(175), 0, pi)]; % 175 degrees
        equality = [equality; false(1, 1)];

        phi_kite = atan2(T_kite_B(3), T_kite_B(2)); % Angle in y-z plane
        g = {g{:}, si_to_scaled(phi_kite, 0, pi)};
        lbg = [lbg; si_to_scaled(deg2rad(5), 0, pi)]; % 5 degrees
        ubg = [ubg; si_to_scaled(deg2rad(175), 0, pi)]; % 175 degrees
        equality = [equality; false(1, 1)];
    end

    % Calculate aerodynamic variables (apparent wind speed, angle of attack, sideslip angle)
    [Va, alpha, beta] = symbolic_funcs.aerodynamics_vars(X0_si, paramsColumnVec);

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

    % Constraint on generated power
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

        % Enforce that F_aero_B and M_aero_B come from the aerodynamic model
        [F_a_B_expr,  M_a_B_expr] = symbolic_funcs.F_M_aero(Xk_si, cs_defl, paramsColumnVec);
        g = {g{:}, (F_aero_B - F_a_B_expr)./ sqrt((u_U(2:4) - u_L(2:4))/2)}; % Force residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];
        g = {g{:}, (M_aero_B - M_a_B_expr) ./ sqrt((u_U(5:7) - u_L(5:7))/2)}; % Moment residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        % Euler integration for ODE
        rhs = vertcat(symbolic_funcs.winch_dynamics(Xk_si, M_winch, Zk_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(Xk_si, Zk_si, F_aero_B, M_aero_B, F_prop_B, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    elseif params_num.ctrl_mode == 2 % Fictitious force and torque control
        M_winch = Uk_si(1);      % Winch moment
        F_fict_B = Uk_si(2:4);   % Fictitious force in body frame (control variable)
        M_fict_B = Uk_si(5:7);   % Fictitious moment in body frame (control variable)
        F_aero_B = zeros(3,1);   % Aerodynamic force is zero in this mode
        M_aero_B = zeros(3,1);   % Aerodynamic moment is zero in this mode

        % Euler integration for ODE with fictitious forces
        rhs = vertcat(symbolic_funcs.winch_dynamics(Xk_si, M_winch, Zk_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(Xk_si, Zk_si, F_fict_B, M_fict_B, 0, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    else % Homotopy control (blend of aerodynamic and fictitious forces)
        M_winch = Uk_si(1);      % Winch moment
        F_fic_a_B = Uk_si(2:4);  % Fictitious/aerodynamic force (depending on Phi)
        M_fic_a_B = Uk_si(5:7);  % Fictitious/aerodynamic moment (depending on Phi)
        cs_defl = Uk_si(8:10);   % Control surface deflections
        F_prop_B = Uk_si(11);    % Propeller force
        F_aero_B = Uk_si(12:14); % Actual aerodynamic force (used for consistency constraints)
        M_aero_B = Uk_si(15:17); % Actual aerodynamic moment (used for consistency constraints)

        % Enforce that F_aero_B and M_aero_B come from the aerodynamic model
        [F_a_B_expr,  M_a_B_expr] = symbolic_funcs.F_M_aero(Xk_si, cs_defl, paramsColumnVec);
        g = {g{:}, (F_aero_B - F_a_B_expr) ./ sqrt((u_U(12:14) - u_L(12:14))/2)}; % Force residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        g = {g{:}, (M_aero_B - M_a_B_expr) ./ sqrt((u_U(15:17) - u_L(15:17))/2)}; % Moment residuals
        lbg = [lbg; zeros(3,1)];
        ubg = [ubg; zeros(3,1)];
        equality = [equality; true(3, 1)];

        % Combine fictitious/aerodynamic forces/moments based on homotopy parameter Phi(1)
        F_M_hat_a_B = Phi(1) .* [F_fic_a_B; M_fic_a_B] + (1 - Phi(1)) .* [F_aero_B; M_aero_B];
        % Apply homotopy to propeller force based on Phi(2)
        F_hat_p_B = Phi(2) .* F_prop_B;

        % Euler integration for ODE with blended forces
        rhs = vertcat(symbolic_funcs.winch_dynamics(Xk_si, M_winch, Zk_si(1), paramsColumnVec), symbolic_funcs.kite_dynamics(Xk_si, Zk_si, F_M_hat_a_B(1:3), F_M_hat_a_B(4:6), F_hat_p_B, paramsColumnVec));
        ode_eq = Xk_si - X0_si - h .* rhs;

    end

    % Constraint on translational acceleration magnitude
    trans_acc = rhs(3:5)' * rhs(3:5); % Squared magnitude of translational acceleration
    g = {g{:}, si_to_scaled(trans_acc, kite_acc_min, (kite_acc_max * params_num.env.g)^2)}; % Scaled constraint
    lbg = [lbg; 0];
    ubg = [ubg; 1];
    equality = [equality; false(1, 1)];

    % Constraint on winch acceleration
    winch_acc = rhs(2); % Winch acceleration from RHS
    winch_acc_scaled = si_to_scaled(winch_acc , winch_acc_min, winch_acc_max); % Scaled winch acceleration
    g = {g{:}, winch_acc_scaled};
    lbg = [lbg; 0];
    ubg = [ubg; 1];

    % Combine algebraic and ordinary differential equation (ODE) constraints
    % These are equality constraints enforcing model dynamics.
    g = {g{:}, vertcat(alg_eq ./sqrt((z_U - z_L)/2), ode_eq ./ sqrt((x_U - x_L))/2)};
    lbg = [lbg; zeros(Nx+Nz,1)];
    ubg = [ubg; zeros(Nx+Nz,1)];
    equality = [equality; true(Nx+Nz, 1)];

    % Calculate and scale angular velocity derivatives and lift force
    omega_dot = si_to_scaled(rhs(6:8), -0.5*ones(3,1), 0.5*ones(3,1)); % Scaled angular accelerations
    lift_force = -F_aero_B(1)*sin(alpha) - F_aero_B(3)*cos(alpha); % Lift force calculation

end