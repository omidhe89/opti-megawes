function ocp_bounds = system_OCP_bounds(params_num)
% SYSTEM_OCP_BOUNDS Defines the lower and upper bounds for the OCP variables.
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
%   This function defines the numerical lower and upper bounds for the states (x),
%   controls (u), algebraic variables (z), and homotopy parameters (phi) used
%   in the Optimal Control Problem (OCP) formulation. It also sets bounds for
%   other important physical quantities like angle of attack, sideslip angle,
%   apparent wind speed, and winch acceleration.
%
%   The bounds are tailored based on the simulation's `path_type` and
%   `ctrl_mode` settings, ensuring the optimization problem adheres to
%   realistic operational limits of the Airborne Wind Energy (AWE) system.
%
% Inputs:
%   params_num - struct, A structure containing numerical parameters of the
%                system, including `sim.path_type`, `winch.radius`,
%                `sim.winding_number`, and `sim.ctrl_mode`.
%
% Outputs:
%   ocp_bounds - struct, A structure containing all defined lower (`_L`) and
%                upper (`_U`) bounds for the OCP variables and other key
%                physical constraints.
%
% See also: system_vars_params, SOLVE_MEGAWES_OCP

% Extract path_type for conditional bounding logic
path_type = params_num.sim.path_type;

%% Known System Limits
% These are fundamental physical or operational limits of the kite system.
tether_length_min = 20;                     % Minimum allowed tether length (m)
tether_length_max = 800;                    % Maximum allowed tether length (m)
tether_speed_min  = -15;                    % Minimum allowed tether speed (m/s), i.e., retraction speed
tether_speed_max  = 15;                     % Maximum allowed tether speed (m/s), i.e., deployment speed
tether_F_max      = 1.664910873257869e+06;  % Maximum allowed tether force (N)

% Control surface deflection limits (radians)
defl_ai_min = -15 * pi/180; % Minimum aileron deflection
defl_ai_max = 15 * pi/180;  % Maximum aileron deflection
defl_el_min = -10 * pi/180; % Minimum elevator deflection
defl_el_max = 10 * pi/180;  % Maximum elevator deflection
defl_ru_min = -10 * pi/180; % Minimum rudder deflection
defl_ru_max = 10 * pi/180;  % Maximum rudder deflection

%% Set NLP Variables Bounds
% Define lower (L) and upper (U) bounds for each NLP variable type.

% --- State Variable Bounds (x) ---
% x1: Winch angle (rad)
x1L = tether_length_min / params_num.winch.radius;
x1U = tether_length_max / params_num.winch.radius;

% x2: Winch angular velocity (rad/s)
x2L = tether_speed_min / params_num.winch.radius;
x2U = tether_speed_max / params_num.winch.radius;

% x3, x4, x5: Kite body-frame velocities (Vel_K_B_x, Vel_K_B_y, Vel_K_B_z) (m/s)
x3L = 35;                   x3U = 100;
x4L = -100;                 x4U = 100;
x5L = -100;                 x5U = 100;

% x6, x7, x8: Kite body-frame angular velocities (Omega_OB_B_p, Omega_OB_B_q, Omega_OB_B_r) (rad/s)
x6L = -50*pi/180;           x6U = 50*pi/180;

if params_num.sim.path_type == 0 % Example: Specific path type might have tighter bounds
    x7L = -40*pi/180;       x7U = 40*pi/180;
    x8L = -40*pi/180;       x8U = 40*pi/180;
else % Default or other path types
    x7L = -50*pi/180;       x7U = 50*pi/180;
    x8L = -50*pi/180;       x8U = 50*pi/180;
end

% x9, x10, x11: Kite body orientation (roll, pitch, yaw angles relative to ground frame) (rad)
x9L  = -75*pi/180;          x9U  = 75*pi/180;  % Roll angle
x10L = -75*pi/180;          x10U = 75*pi/180;  % Pitch angle

if path_type == 8 % Specific path type, e.g., lemniscates
    x11L = -pi + pi/6;      x11U = pi - pi/6; % Yaw angle for lemniscate path
else % Default or circular paths (CW/CCW)
    
    % Initialize bounds for circular paths (CW)
    x11U = (2 * params_num.sim.winding_number + 0.05) * pi + pi/2;
    x11L = -1.05 * pi/2; 
    % Re-assign bounds for circular paths (CCW) by symmetry
    x11L = -x11U;
    x11U = -x11L;
end

% x12, x13, x14: Aircarft position in spherical coordinates (longitude, latitude, height)
x12L = -pi/2 + pi/20;       x12U = pi/2 - pi/20;   % Longitude (azimuth) relative to ground (rad)
x13L = pi/20;               x13U = pi/2 - (pi/20); % Latitude (elevation) relative to ground (rad)
x14L = 50;                  x14U = 800;            % Kite altitude (h_tau) (m) - bounds should be verified for consistency with tether length

% Aggregate state bounds into vectors
x_L = [x1L; x2L; x3L; x4L; x5L; x6L; x7L; x8L; x9L; x10L; x11L; x12L; x13L; x14L];
x_U = [x1U; x2U; x3U; x4U; x5U; x6U; x7U; x8U; x9U; x10U; x11U; x12U; x13U; x14U];

% --- Control Variable Bounds (u) ---
% The structure and number of control variables depend on the control mode.
if params_num.ctrl_mode == 1
    % Control Mode 1: Winch torque and control surfaces deflection are physical control inputs.
    %                 The aerodynamic forces and torques are included as additional terms to reduce symbolic complexity.
    u1L = 1.31e3;                             u1U = 2.1 * params_num.winch.radius * tether_F_max; % Winch torque command (Nm)
    u2L = -10e5;                              u2U = 10e5;   % F_a_B(1) - Aerodynamic force along body X-axis (N)
    u3L = -10e5;                              u3U = 10e5;   % F_a_B(2) - Aerodynamic force along body Y-axis (N)
    u4L = -30e5;                              u4U = 30e5;   % F_a_B(3) - Aerodynamic force along body Z-axis (N)
    u5L = -26e5;                              u5U = 26e5;   % M_a_B(1) - Aerodynamic moment about body X-axis (Nm)
    u6L = -16e5;                              u6U = 16e5;   % M_a_B(2) - Aerodynamic moment about body Y-axis (Nm)
    u7L = -16e5;                              u7U = 16e5;   % M_a_B(3) - Aerodynamic moment about body Z-axis (Nm)
    u8L = defl_ai_min;                        u8U = defl_ai_max;  % Aileron deflection (rad)
    u9L = defl_el_min;                        u9U = defl_el_max;  % Elevator deflection (rad)
    u10L = defl_ru_min;                       u10U = defl_ru_max; % Rudder deflection (rad)
    u11L = 0;                                 u11U = 0.5e-2;      % Propulsive force (N) almost zero let kite to glide.
    u_L = [u1L; u2L; u3L; u4L; u5L; u6L; u7L; u8L; u9L; u10L; u11L];
    u_U = [u1U; u2U; u3U; u4U; u5U; u6U; u7U; u8U; u9U; u10U; u11U];
elseif params_num.ctrl_mode == 2
    % Control Mode 2: Winch torque and fictitious forces/torques
    u1L = 1.31e3;                             u1U = 2.1 * params_num.winch.radius * tether_F_max; % Winch torque command (Nm)
    u2L = -5e5;                               u2U = 5e5;    % F_fict_B(1) - Fictitious force along body X-axis (N)
    u3L = -5e5;                               u3U = 5e5;    % F_fict_B(2) - Fictitious force along body Y-axis (N)
    u4L = -10e5;                              u4U = 10e5;   % F_fict_B(3) - Fictitious force along body Z-axis (N)
    u5L = -3e5;                               u5U = 3e5;    % M_fict_B(1) - Fictitious moment about body X-axis (Nm)
    u6L = -9e4;                               u6U = 9e4;    % M_fict_B(2) - Fictitious moment about body Y-axis (Nm)
    u7L = -9e4;                               u7U = 9e4;    % M_fict_B(3) - Fictitious moment about body Z-axis (Nm)
    u_L = [u1L; u2L; u3L; u4L; u5L; u6L; u7L];
    u_U = [u1U; u2U; u3U; u4U; u5U; u6U; u7U];
else % Control Mode 3: Winch torque and a combination of deflections and fictitious forces
     %                 The aerodynamic forces and torques are included as additional terms to reduce symbolic complexity.
    u1L = 1.31e3;                             u1U = 2.1 * params_num.winch.radius * tether_F_max; % Winch torque command (Nm)
    u2L = -2e5;                               u2U = 2e5;    % F_fict_B(1) - Fictitious force along body X-axis (N)
    u3L = -2e5;                               u3U = 2e5;    % F_fict_B(2) - Fictitious force along body Y-axis (N)
    u4L = -6e5;                               u4U = 6e5;    % F_fict_B(3) - Fictitious force along body Z-axis (N)
    u5L = -1.5e5;                             u5U = 1.5e5;  % M_fict_B(1) - Fictitious moment about body X-axis (Nm)
    u6L = -6e4;                               u6U = 6e4;    % M_fict_B(2) - Fictitious moment about body Y-axis (Nm)
    u7L = -6e4;                               u7U = 6e4;    % M_fict_B(3) - Fictitious moment about body Z-axis (Nm)
    u8L = defl_ai_min;                        u8U = defl_ai_max;  % Aileron deflection (rad)
    u9L = defl_el_min;                        u9U = defl_el_max;  % Elevator deflection (rad)
    u10L = defl_ru_min;                       u10U = defl_ru_max; % Rudder deflection (rad)
    u11L = -0.5e3;                            u11U = 0.5e3;       % Propulsive force (N)
    u12L = -10e5;                             u12U = 10e5;  % F_a_B(1) - Aerodynamic force along body X-axis (N)
    u13L = -10e5;                             u13U = 10e5;  % F_a_B(2) - Aerodynamic force along body Y-axis (N)
    u14L = -40e5;                             u14U = 40e5;  % F_a_B(3) - Aerodynamic force along body Z-axis (N)
    u15L = -26e5;                             u15U = 26e5;  % M_a_B(1) - Aerodynamic moment about body X-axis (Nm)
    u16L = -16e5;                             u16U = 16e5;  % M_a_B(2) - Aerodynamic moment about body Y-axis (Nm)
    u17L = -16e5;                             u17U = 16e5;  % M_a_B(3) - Aerodynamic moment about body Z-axis (Nm)
    u_L = [u1L; u2L; u3L; u4L; u5L; u6L; u7L; u8L; u9L; u10L; u11L; u12L; u13L; u14L; u15L; u16L; u17L];
    u_U = [u1U; u2U; u3U; u4U; u5U; u6U; u7U; u8U; u9U; u10U; u11U; u12U; u13U; u14U; u15U; u16U; u17U];
end

% --- Algebraic Variable Bounds (z) ---
% These represent auxiliary variables, often forces or geometric properties.
z1L = 1.31e3; z1U = tether_F_max;              % Tether force (N)
z2L = pi/50;  z2U = pi/2 - (pi/50);            % Example: Tether azimuth angle (rad)
z3L = -pi/2;  z3U = pi/2;                      % Example: Tether elevation angle (rad)

% Aggregate algebraic bounds into vectors
z_L = [z1L; z2L; z3L];
z_U = [z1U; z2U; z3U];

% --- Homotopy Parameter Bounds (phi) ---
% These parameters are typically used to smoothly transition between different
% problem formulations or objective functions during optimization.
phi1L = 1; phi1U = 1; % Controls switch from fictitious to aerodynamic forces.
phi2L = 1; phi2U = 1; % Controls reduction of Propulsive force and initiation of gliding.
phi3L = 1; phi3U = 1; % Controls maximization of power.


% Aggregate homotopy parameter bounds into vectors
phi_L = [phi1L; phi2L; phi3L];
phi_U = [phi1U; phi2U; phi3U];

% --- Time Bounds ---
% Defines the allowable range for the optimal control problem's time horizon.
time_L = params_num.sim.Tf - 0.75 * (params_num.sim.Tf / params_num.sim.winding_number);
time_U = params_num.sim.Tf + 0.75 * (params_num.sim.Tf / params_num.sim.winding_number);

%% Additional constraints
% aerodynamaics
alpha_min     = -15 * pi/180;       alpha_max     = 5 * pi/180;   % Angle of attack (rad)
beta_min      = -10 * pi/180;       beta_max      = 10 * pi/180;  % Sideslip angle (rad)
v_app_min     = 20;                 v_app_max     = 120;          % Apparent wind speed (m/s)
% accelerations 
winch_acc_min = -1.6;               winch_acc_max = 1.6;          % Winch acceleration (m/s^2)
kite_acc_min = 0;                   kite_acc_max = 4;             % kite linear acceleration magnitute ( times g=9.81 m/s^2)

% flight region
pos_W_x_min = 0;                    pos_W_x_max = 1000;
if params_num.sim.path_type == 0
    pos_W_y_min = -250;             pos_W_y_max = 250;  % circular trajectory 
else
    pos_W_y_min = -400;             pos_W_y_max = 400; % lemniscate trajectory
end
pos_W_z_min = 150;                  pos_W_z_max = 600;

% power 
power_min = -3e6;                   power_max = 3e6; % harvested power
%% Populate Output Structure
% Consolidate all defined bounds into a single output structure.
ocp_bounds = struct;
ocp_bounds.x_L = x_L;
ocp_bounds.x_U = x_U;
ocp_bounds.u_L = u_L;
ocp_bounds.u_U = u_U;
ocp_bounds.z_L = z_L;
ocp_bounds.z_U = z_U;
ocp_bounds.phi_L = phi_L;
ocp_bounds.phi_U = phi_U;
ocp_bounds.alpha_min = alpha_min;
ocp_bounds.alpha_max = alpha_max;
ocp_bounds.beta_min = beta_min;
ocp_bounds.beta_max = beta_max;
ocp_bounds.v_app_min = v_app_min;
ocp_bounds.v_app_max = v_app_max;
ocp_bounds.winch_acc_min = winch_acc_min;
ocp_bounds.winch_acc_max = winch_acc_max;
ocp_bounds.time_L = time_L;
ocp_bounds.time_U = time_U;
ocp_bounds.kite_acc_min = kite_acc_min;
ocp_bounds.kite_acc_max = kite_acc_max;
ocp_bounds.pos_W_x_min = pos_W_x_min;
ocp_bounds.pos_W_x_max = pos_W_x_max;
ocp_bounds.pos_W_y_min = pos_W_y_min;
ocp_bounds.pos_W_y_max = pos_W_y_max;
ocp_bounds.pos_W_z_min = pos_W_z_min;
ocp_bounds.pos_W_z_max = pos_W_z_max;
ocp_bounds.power_min = power_min;
ocp_bounds.power_max = power_max;

end % End of function system_OCP_bounds