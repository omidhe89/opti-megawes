function [vars, params_sym, params_num] = system_vars_params(windSpeed, windDir, params_num)
% SYSTEM_VARS_PARAMS Defines symbolic variables and system parameters for MegAWES OCP.
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
%   This function initializes all symbolic variables (states, controls,
%   algebraic variables, derivatives, and homotopy parameters) using CasADi's
%   SX.sym. It also defines both the symbolic (`params_sym`) and numerical
%   (`params_num`) values for various system parameters, such as environment
%   properties, winch characteristics, tether properties, and kite
%   aerodynamics. The structure of control inputs (`Nu`) varies based on the
%   selected control mode (`ctrl_mode`), and aerodynamic coefficients are
%   loaded or set based on the `aero_model` choice.
%
% Inputs:
%   windSpeed  - double, Base wind speed at a reference height (m/s).
%   windDir    - double, Base wind direction (radians).
%   params_num - struct, Numerical parameters structure, partially populated
%                from the main script, e.g., containing `sim.ctrl_mode` and
%                `sim.aero_model`, `homotopy.phi_numbers`.
%
% Outputs:
%   vars       - struct, Contains all defined symbolic variables (`SX.sym` objects)
%                for states (x), controls (u), algebraic variables (z),
%                derivatives (x_dot), homotopy parameters (phi), and various
%                force/torque/deflection variables. It also includes symbolic
%                representations of variable bounds.
%   params_sym - struct, Contains all defined symbolic parameters (`SX.sym` objects)
%                for environment, winch, tether, and kite properties.
%   params_num - struct, Updated numerical parameters structure, filling with
%                default or loaded numerical values corresponding to the
%                symbolic parameters.
%
% Dependencies:
%   - CasADi library (imported)
%   - DE2019_params.mat (data file containing aircraft parameters, copied from
%     the MegAWES repository under the Apache-2.0 license)
%   - set_borne_aero_params_num.m (function to set aerodynamic coefficients)
%
% See also: SOLVE_MEGAWES_OCP, build_dynamics

% Import CasADi library for symbolic expression generation.
import casadi.*;
load DE2019_params.mat DE2019 % Load numerical parameters for the kite (aircraft) from a MAT file.

%% --- Define System Dimensions (Number of States, Controls, Algebraic Variables) ---

% Nx: Number of differential states.
% x(1:2) typically represents winch-related states (e.g., roration angle, rotation speed).
% x(3:14) typically represent aircraft states (e.g., velocity, angular rates, orientation, position(in spherical refernece frame)).
params_num.Nx = 14;  

% Nu: Number of control inputs. This varies significantly based on the control mode.
if params_num.ctrl_mode == 1 % Aerodynamic Control Mode
    % u(1): Winch control input.
    % u(2:4): Aerodynamic forces on the kite in body frame (F_a_B).
    % u(5:7): Aerodynamic moments on the kite in body frame (M_a_B).
    % u(8:10): Control surface deflections (e.g., aileron, elevator, rudder).
    % u(11): Propulsive force on the kite in body frame (F_p_B).
    params_num.Nu = 11;
elseif params_num.ctrl_mode == 2 % Fictitious Control Mode
    % u(1): Winch control input.
    % u(2:4): Fictitious forces on the kite in body frame (F_fict_B).
    % u(5:7): Fictitious moments on the kite in body frame (M_fict_B).
    % Fictitious controls are often used for simpler problem formulations or
    % as intermediate steps in homotopy methods.
    params_num.Nu = 7;
else % Homotopy Control Mode (params_num.ctrl_mode == 3)
    % This mode uses a combination of fictitious controls, control surface
    % deflections, propulsive force, and possibly direct aerodynamic forces/moments.
    % The high number of controls suggests a rich set of control variables
    % that are gradually introduced or constrained by the homotopy method.
    % u(1): Winch control input.
    % u(2:4): Fictitious forces on the kite in body frame (F_fict_B).
    % u(5:7): Fictitious moments on the kite in body frame (M_fict_B).
    % u(8:10): Control surface deflections.
    % u(11): Propulsive force.
    % u(12:14): Aerodynamic forces on the kite in body frame (F_aero_B).
    % u(15:17): Aerodynamic moments on the kite in body frame (M_aero_B).
    params_num.Nu = 17;
end

% Nz: Number of algebraic states.
% These are states whose dynamics are described by algebraic equations.
% z(1:3): Direction (elevation, azimuth) and magnitude of tether force at the winch.
params_num.Nz = 3;

%% --- Initialize Symbolic Variables (`vars` Structure) ---
% Creates CasADi symbolic variables (SX.sym) that represent the dynamic
% states, control inputs, and algebraic variables throughout the OCP.
% These variables will be the decision variables for the NLP solver.
vars = struct();

% `x`: Differential states vector.
vars.x = SX.sym('x', params_num.Nx);
% `u`: Control inputs vector.
vars.u = SX.sym('u', params_num.Nu);
% `z`: Algebraic states vector.
vars.z = SX.sym('z', params_num.Nz);
% `x_dot`: Derivatives of the differential states (not used in dynamics).
vars.x_dot = SX.sym('x_dot', params_num.Nx);
% `phi`: Homotopy parameters. These are typically scalar parameters that
% linearly blend between a simpler problem and the full complex problem.
vars.phi = SX.sym('phi', params_num.homotopy.phi_numbers);


% Individual symbolic variables for specific forces, torques, and deflections.
% These might be used internally by dynamics functions or as components of `vars.u`.
vars.kite_aero_force = SX.sym('F_a_B', 3);      % Aerodynamic force on kite in body frame.
vars.kite_prop_force = SX.sym('F_p_B', 1);      % Propulsive force on kite in body frame.
vars.kite_aero_torque = SX.sym('M_a_B', 3);     % Aerodynamic torque on kite in body frame.
vars.kite_cs_deflections = SX.sym('defl', 3);   % Control surface deflections.
vars.winch_torque = SX.sym('M_w', 1);           % Winch torque.

%% --- Initialize Symbolic and Numerical Parameters (`params_sym`, `params_num`) ---
% Defines both symbolic placeholders and their corresponding numerical values
% for physical and environmental parameters of the system.
params_sym = struct();

% Environment Parameters
params_sym.env.g = SX.sym('g');               % Gravitational acceleration (symbolic).
params_num.env.g = 9.81;                      % Gravitational acceleration (numerical, m/s^2).
params_sym.env.wind.direction = SX.sym('wind_dir'); % Wind direction (symbolic).
params_num.env.wind.direction = windDir;      % Wind direction (numerical, radians), passed from input.
params_sym.env.wind.base_velocity = SX.sym('wind_velocity'); % Base wind velocity (symbolic).
params_num.env.wind.base_velocity = windSpeed;% Base wind velocity (numerical, m/s), passed from input.
params_sym.env.wind.rho = SX.sym('wind_rho'); % Air density (symbolic).
params_num.env.wind.rho = 1.225;              % Air density (numerical, kg/m^3).
params_sym.env.wind.z0 = SX.sym('wind_z0');   % Surface roughness length (symbolic).
params_num.env.wind.z0 = 0.002;               % Surface roughness length (numerical, m).
params_sym.env.wind.z_ref = SX.sym('wind_z_ref'); % Reference height for wind speed (symbolic).
params_num.env.wind.z_ref = 100;              % Reference height for wind speed (numerical, m).

% Winch Parameters
params_sym.winch.radius = SX.sym('winch_radius'); % Winch drum radius (symbolic).
params_num.winch.radius = 1.5;                    % Winch drum radius (numerical, m).
params_sym.winch.inertia = SX.sym('winch_inertia'); % Winch system inertia (symbolic).
params_num.winch.inertia = 4033 + 1329 + 1786;   % Winch + flange + rotor inertia (numerical, kg*m^2).
params_sym.winch.friction = SX.sym('winch_friction'); % Winch friction coefficient (symbolic).
params_num.winch.friction = 10;                   % Winch friction coefficient (numerical).

% Tether Parameters
params_sym.tether.Np = SX.sym('tether_Np');   % Number of tether segments (symbolic).
params_num.tether.Np = 11;                    % Number of tether segments (numerical).
params_sym.tether.rho = SX.sym('tether_rho'); % Tether material density (symbolic).
params_num.tether.rho = 0.672901477941722;    % Tether material density (numerical, kg/m).
params_sym.tether.E = SX.sym('tether_E');     % Tether Young's Modulus (elasticity) (symbolic).
params_num.tether.E = 1.16e+11;               % Tether Young's Modulus (numerical, Pa).
params_sym.tether.A = SX.sym('tether_A');     % Tether cross-sectional area (symbolic).
params_num.tether.A = 6.937128638574452e-04;  % Tether cross-sectional area (numerical, m^2).
params_sym.tether.d = SX.sym('tether_d');     % Tether diameter (symbolic).
params_num.tether.d = 0.029719735041800;      % Tether diameter (numerical, m).
params_sym.tether.cd = SX.sym('tether_cd');   % Tether drag coefficient (symbolic).
params_num.tether.cd = 1.1;                   % Tether drag coefficient (numerical).

% Kite (Aircraft) Parameters
params_sym.kite.mass = SX.sym('kite_mass');   % Kite mass (symbolic).
params_num.kite.mass = DE2019.mass;           % Kite mass (numerical, kg).
params_sym.kite.Jx = SX.sym('kite_Jx');       % Kite moment of inertia around x-axis (symbolic).
params_num.kite.Jx = DE2019.Jx;               % Kite moment of inertia around x-axis (numerical, kg*m^2).
params_sym.kite.Jz = SX.sym('kite_Jz');       % Kite moment of inertia around z-axis (symbolic).
params_num.kite.Jz = DE2019.Jz;               % Kite moment of inertia around z-axis (numerical, kg*m^2).
params_sym.kite.Jxz = SX.sym('kite_Jxz');     % Kite product of inertia Jxz (symbolic).
params_num.kite.Jxz = DE2019.Jxz;             % Kite product of inertia Jxz (numerical, kg*m^2).
params_sym.kite.Jy = SX.sym('kite_Jy');       % Kite moment of inertia around y-axis (symbolic).
params_num.kite.Jy = DE2019.Jy;               % Kite moment of inertia around y-axis (numerical, kg*m^2).

% Wing Parameters (part of Kite)
params_sym.kite.b = SX.sym('kite_b');         % Wing span (symbolic).
params_num.kite.b = DE2019.b;                 % Wing span (numerical, m).
params_sym.kite.c = SX.sym('kite_c');         % Mean aerodynamic chord (symbolic).
params_num.kite.c = DE2019.c;                 % Mean aerodynamic chord (numerical, m).
params_sym.kite.s_ref = SX.sym('kite_S_ref'); % Reference wing area (symbolic).
params_num.kite.s_ref = DE2019.S_ref;         % Reference wing area (numerical, m^2).


% Aerodynamic Model Specific Parameters
if strcmp(params_num.sim.aero_model, 'MegAWES')
    % These parameters are specific to the 'MegAWES' aerodynamic model,
    % which involves a detailed geometric and mass distribution model.
    % More info: 
    %       D. Eijkelhof, R. Schmehl: Six-degrees-of-freedom simulation model for future multi-megawatt airborne wind energy systems.
    %       Renewable Energy, Vol. 196, pp. 137-150, 2022. doi:10.1016/j.renene.2022.06.094

    params_sym.kite.cog = SX.sym('kite_cog', 3); % Center of gravity coordinates (symbolic).
    params_num.kite.cog = DE2019.cog;            % Center of gravity coordinates (numerical, m).
    params_sym.kite.wing_cRoll_Static = SX.sym('wing_cRoll_Static'); % Static roll coefficient of the wing (symbolic).
    params_num.kite.wing_cRoll_Static = DE2019.initAircraft.wing_cRoll_Static; % Static roll coefficient (numerical).

    % Elevator (Horizontal Stabilizer) Parameters
    params_sym.kite.s_el = SX.sym('kite_S_el');   % Elevator reference area (symbolic).
    params_num.kite.s_el = DE2019.S_el;           % Elevator reference area (numerical, m^2).
    params_sym.kite.c_el = SX.sym('kite_c_el');   % Elevator chord (symbolic).
    params_num.kite.c_el = DE2019.c_el;           % Elevator chord (numerical, m).
    params_sym.kite.b_el = SX.sym('kite_b_el');   % Elevator span (symbolic).
    params_num.kite.b_el = DE2019.b_el;           % Elevator span (numerical, m).
    params_sym.kite.x_el = SX.sym('kite_x_el');   % Elevator x-position relative to reference (symbolic).
    params_num.kite.x_el = DE2019.x_el;           % Elevator x-position (numerical, m).

    % Rudder (Vertical Stabilizer) Parameters
    params_sym.kite.s_vs = SX.sym('kite_S_vs');   % Rudder reference area (symbolic).
    params_num.kite.s_vs = DE2019.S_vs;           % Rudder reference area (numerical, m^2).
    params_sym.kite.c_vs = SX.sym('kite_c_vs');   % Rudder chord (symbolic).
    params_num.kite.c_vs = DE2019.c_vs;           % Rudder chord (numerical, m).
    params_sym.kite.b_vs = SX.sym('kite_b_vs');   % Rudder span (symbolic).
    params_num.kite.b_vs = DE2019.b_vs;           % Rudder span (numerical, m).
    params_sym.kite.x_vs = SX.sym('kite_x_vs');   % Rudder x-position relative to reference (symbolic).
    params_num.kite.x_vs = DE2019.x_vs;           % Rudder x-position (numerical, m).
else
    % If not 'MegAWES' model, assume a one of aerodynamic models (ALM, VLM, CFD) calculated
    % in BORNE project is used.
  
    % These coefficients are typically functions of angle of attack (alpha),
    % sideslip angle (beta), and control surface deflections.

    % CX Components
    params_sym.kite.CX = struct();
    params_sym.kite.CX.zero = SX.sym('CX0');          % Zero-lift drag coefficient.
    params_sym.kite.CX.alpha = SX.sym('CX_alpha', 2); % Coefficients for alpha dependency.
    params_sym.kite.CX.q = SX.sym('CX_q', 3);         % Coefficients for pitch rate (q) dependency.
    params_sym.kite.CX.delta_e = SX.sym('CX_delta_e', 3); % Coefficients for elevator deflection.

    % CY Components
    params_sym.kite.CY = struct();
    params_sym.kite.CY.beta = SX.sym('CY_beta', 3);   % Coefficients for sideslip angle (beta) dependency.
    params_sym.kite.CY.p = SX.sym('CY_p', 3);         % Coefficients for roll rate (p) dependency.
    params_sym.kite.CY.r = SX.sym('CY_r', 3);         % Coefficients for yaw rate (r) dependency.
    params_sym.kite.CY.delta_a = SX.sym('CY_delta_a', 3); % Coefficients for aileron deflection.
    params_sym.kite.CY.delta_r = SX.sym('CY_delta_r', 3); % Coefficients for rudder deflection.

    % CZ Components
    params_sym.kite.CZ = struct();
    params_sym.kite.CZ.zero = SX.sym('CZ0');         % Zero-lift coefficient.
    params_sym.kite.CZ.alpha = SX.sym('CZ_alpha', 2); % Coefficients for alpha dependency.
    params_sym.kite.CZ.q = SX.sym('CZ_q', 3);         % Coefficients for pitch rate (q) dependency.
    params_sym.kite.CZ.delta_e = SX.sym('CZ_delta_e', 3); % Coefficients for elevator deflection.

    % Cl (Rolling Moment Coefficient) Components
    params_sym.kite.Cl = struct();
    params_sym.kite.Cl.beta = SX.sym('Cl_beta', 3);   % Coefficients for beta dependency.
    params_sym.kite.Cl.p = SX.sym('Cl_p', 3);         % Coefficients for roll rate (p) dependency.
    params_sym.kite.Cl.r = SX.sym('Cl_r', 3);         % Coefficients for yaw rate (r) dependency.
    params_sym.kite.Cl.delta_a = SX.sym('Cl_delta_a', 3); % Coefficients for aileron deflection.
    params_sym.kite.Cl.delta_r = SX.sym('Cl_delta_r', 3); % Coefficients for rudder deflection.

    % Cm (Pitching Moment Coefficient) Components
    params_sym.kite.Cm = struct();
    params_sym.kite.Cm.zero = SX.sym('Cm0');         % Zero-pitching moment coefficient.
    params_sym.kite.Cm.alpha = SX.sym('Cm_alpha', 2); % Coefficients for alpha dependency.
    params_sym.kite.Cm.q = SX.sym('Cm_q', 3);         % Coefficients for pitch rate (q) dependency.
    params_sym.kite.Cm.delta_e = SX.sym('Cm_delta_e', 3); % Coefficients for elevator deflection.

    % Cn (Yawing Moment Coefficient) Components
    params_sym.kite.Cn = struct();
    params_sym.kite.Cn.beta = SX.sym('Cn_beta', 3);   % Coefficients for beta dependency.
    params_sym.kite.Cn.p = SX.sym('Cn_p', 3);         % Coefficients for roll rate (p) dependency.
    params_sym.kite.Cn.r = SX.sym('Cn_r', 3);         % Coefficients for yaw rate (r) dependency.
    params_sym.kite.Cn.delta_a = SX.sym('Cn_delta_a', 3); % Coefficients for aileron deflection.
    params_sym.kite.Cn.delta_r = SX.sym('Cn_delta_r', 3); % Coefficients for rudder deflection.

    % Call a separate function to populate the numerical values for the BORNE
    % aerodynamic coefficients based on the selected aerodynamic model.
    [params_num.kite.CX, params_num.kite.CY, params_num.kite.CZ, ...
     params_num.kite.Cl, params_num.kite.Cm, params_num.kite.Cn] = set_borne_aero_params_num(params_num.sim.aero_model);
end

end % End of function system_vars_params

