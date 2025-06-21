function [kite_data, tether_dynamics, tether_alg_equation, kite_dynamics, winch_dynamics, aerodymaics_vars, F_M_aero, params_array] = build_dynamics(vars, params_sym, params_num)
% BUILD_DYNAMICS Generates symbolic expressions for an Airborne Wind Energy System (AWES).
%
% Copyright (c) 2025 Omid Heydarnia
%
% This function is part of the 'Opti-MegAWES' tool and is licensed under the
% Apache License, Version 2.0. See the accompanying LICENSE file for the
% full text of the license.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 12, 2025
% Last Modified: June 12, 2025
%
% Description:
%   This function constructs symbolic representations of the dynamics for an
%   Airborne Wind Energy System (AWES) using CasADi. It defines the relationships
%   between states (x), controls (u), algebraic variables (z), and system
%   parameters (p) to derive the time derivatives (x_dot) and algebraic
%   constraints.
%
%   It integrates various sub-models including:
%   - Kite kinematic and dynamic calculations.
%   - Quasi-static tether model.
%   - Winch dynamics.
%   - Aerodynamic force and moment generation based on selected aerodynamic models.
%
% Inputs:
%   vars       - struct, Contains symbolic CasADi variables for states (x),
%                controls (u), algebraic variables (z), and specific control
%                components (F_a_B, M_a_B, F_p_B, deflections, mechanical_winch_torque).
%   params_sym - struct, Contains symbolic system parameters. 
%   params_num - struct, Contains numerical simulation parameters.
%
% Outputs:
%   kite_data           - CasADi Function for kite kinematic data e.g. kite position and velocities in wind reference frame.
%   tether_dynamics     - CasADi Function for tether model outputs (algebraic function, tether tension, mass positions).
%   tether_alg_equation - CasADi Function for the tether's algebraic constraint.
%   kite_dynamics       - CasADi Function for the kite's state derivatives (dx/dt).
%   winch_dynamics      - CasADi Function for the winch's state derivatives (dx/dt).
%   aerodymaics_vars    - CasADi Function for derived aerodynamic variables (Va, alpha, beta).
%   F_M_aero            - CasADi Function for aerodynamic forces and moments (F_a_B, M_a_B).
%   params_array        - Cell array of symbolic parameters, useful for CasADi function generation.
%
% See also: estimations_ap4, tether_quasi_static_single_shooting,
%           tether_quasi_static_single_shooting_rigid, calc_ap4_dynamic,
%           calc_winch_dynamics, wind_shear, transformFromWtoO,
%           aero_MegAWES, aero_borne, struct_to_key_value_array

    import casadi.*; % Import CasADi functionalities

    % Assign symbolic variables from the input 'vars' structure for clarity
    x = vars.x; % States vector
    u = vars.u; % Controls vector (may contain various components depending on ctrl_mode)
    z = vars.z; % Algebraic variables vector

    % Extract specific control/force components for clearer use
    F_a_B               = vars.kite_aero_force;         % Aerodynamic force in body frame
    M_a_B               = vars.kite_aero_torque;        % Aerodynamic moment in body frame
    F_p_B               = vars.kite_prop_force;         % Propulsive/additional force in body frame (e.g., from propellor)
    deflections         = vars.kite_cs_deflections;     % Control surface deflections (aileron, elevator, rudder)
    mechanical_winch_torque = vars.winch_torque;        % Mechanical torque applied to the winch drum

    % Convert the symbolic parameters structure to a key-value cell array
    % This format is required by CasADi functions.
    params_array = struct_to_key_value_array(params_sym);
    % Extract keys and values separately for CasADi Function inputs
    params_array_keys = vertcat(params_array(:,1));
    params_array_values = vertcat(params_array{:,2});

    % --- Kite Kinematics and Position Estimations ---
    % Calculate various kite-related position and velocity vectors.
    % 'pos_O': Kite position in ground-fixed frame.
    % 'pos_W': Kite position relative to winch in wind frame.
    % 'v_k_W': Kite velocity relative to winch in wind frame.
    % 'M_OB': Rotation matrix from body frame to ground frame.
    % 'v_k_tau': Kite velocity in tangentional refernce frame (for more info refer to Sebastian Rapp's theis).
    [pos_O, pos_W, v_k_W, M_OB, v_k_tau] = estimations_ap4(x, params_sym);
    % Create a CasADi Function for these kite data outputs
    kite_data = Function('kite_data', {x, params_array_values}, {pos_O, pos_W, v_k_W, M_OB, v_k_tau}, ...
                         {'x', 'p'}, {'pos_O', 'pos_W', 'v_k_W', 'M_OB', 'v_k_tau'});

    % --- Tether Dynamics ---
    % Calculate tether tension and other related quantities.
    tether_length = params_sym.winch.radius * x(1); % Tether length derived from winch angle

    if params_num.sim.ablation_study
        % Use a rigid tether model for ablation studies (simplified tether dynamics).
        [alg_fun, T_kite, mass_pos] = tether_quasi_static_single_shooting_rigid(z, [pos_W; v_k_W], tether_length, params_sym, params_num.tether.Np);
    else
        % Use the full quasi-static tether model.
        [alg_fun, T_kite, mass_pos] = tether_quasi_static_single_shooting(z, [pos_W; v_k_W], tether_length, params_sym, params_num.tether.Np);
    end
    % CasADi Function for all tether dynamics outputs
    tether_dynamics = Function('tether_dyn', {x, z, params_array_values}, {alg_fun, T_kite, mass_pos}, ...
                               {'x', 'z', 'p'}, {'alg_fun', 'T_kite', 'mass_pos'});
    % CasADi Function specifically for the algebraic equation (constraint) from the tether model
    tether_alg_equation = Function('tether_alg_equ', {x, z, params_array_values}, {alg_fun}, ...
                                   {'x', 'z', 'p'}, {'alg_fun'});

    % --- Kite Body Dynamics ---
    % Calculate the time derivatives of the kite's states (e.g., velocities, angular velocities, orientations).
    kite_dx_expr = calc_ap4_dynamic(x, [F_a_B, M_a_B], F_p_B, v_k_W, T_kite, params_sym);
    % Create a CasADi Function for the kite's dynamics (state derivatives)
    kite_dynamics = Function('kite_dyn', {x, z, F_a_B, M_a_B, F_p_B, params_array_values}, {kite_dx_expr}, ...
                             {'x', 'z', 'F_a_B', 'M_a_B', 'F_p_B', 'p'}, {'kite_xdot_expr'});

    % --- Winch Dynamics ---
    % Calculate the time derivatives of the winch's states (winch angle, angular velocity).
    tether_force_at_winch = z(1); % Tether force at the winch
    winch_dx_expr = calc_winch_dynamics(x, tether_force_at_winch, mechanical_winch_torque, params_sym);
    % Create a CasADi Function for the winch's dynamics (state derivatives)
    winch_dynamics = Function('winch_dyn', {x, mechanical_winch_torque, tether_force_at_winch, params_array_values}, {winch_dx_expr}, ...
                              {'x', 'mechanical_winch_torque', 'tether_force_at_winch', 'p'}, {'winch_xdot_expr'});

    % --- Aerodynamics Variables and Model Selection ---
    % Calculate aerodynamic relevant variables such as apparent wind speed, angle of attack, and sideslip angle.
    altitude = -pos_O(3); % Altitude of the kite (assuming pos_O(3) is Z-coordinate in ground frame)

    % Calculate wind velocity at kite's position
    v_w_W = wind_shear(M_OB, altitude, params_sym); % Wind velocity in world frame
    v_a_W = v_k_W - v_w_W; % Apparent wind velocity in world frame
    Va    = norm(v_a_W);   % Apparent wind speed

    % Transform apparent wind velocity to kite body frame for angle calculations
    v_a_O = transformFromWtoO(params_sym.env.wind.direction, v_a_W); 
    v_a_B = M_OB' * v_a_O; 

    % Calculate angle of attack (alpha) and sideslip angle (beta)
    alpha = (atan2( v_a_B(3), v_a_B(1) )); % Angle of attack
    beta  = (atan2( v_a_B(2), v_a_B(1) )); % Sideslip angle

    % Air density factor (currently fixed to 1, could incorporate altitude effects)
    % H_rho = 8550; % in meters (8.55 km) - typical scale height for air density
    % rho_factor = exp(altitude / (2 * H_rho)); % Air density variation with altitude
    rho_factor = 1; % For now, assume constant air density

    % Create a CasADi Function for the calculated aerodynamic variables
    aerodymaics_vars = Function('aero_variables', {x, params_array_values}, {Va, alpha, beta}, ...
                                {'x', 'p'}, {'Va', 'alpha', 'beta'});

    % Select and call the appropriate aerodynamic force and moment model
    if strcmp(params_num.sim.aero_model, 'MegAWES')
        % Use the 'MegAWES' detailed aerodynamic model
        [F_a_B_expr, M_a_B_expr, AeroAngels_deg] = aero_MegAWES(x(6:8), v_a_W, alpha, beta, rho_factor, deflections, params_sym);
    else
        % Use the 'BORNE' (or simplified) aerodynamic model based on coefficients
        [F_a_B_expr, M_a_B_expr, AeroAngels_deg] = aero_borne(x(6:8), v_a_W, alpha, beta, rho_factor, deflections, params_sym);
    end

    % Create a CasADi Function for the aerodynamic forces and moments,
    % which depend on states (x), control deflections (u), and parameters (p).
    F_M_aero = Function('aero_model', {x, deflections, params_array_values}, {F_a_B_expr, M_a_B_expr}, ...
                        {'x', 'u', 'p'}, {'F_aero_B', 'M_aero_B'});

end % End of function build_dynamics