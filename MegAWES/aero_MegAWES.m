% Copyright (c) 2025 Omid Heydarnia
%
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
%
% Author: Omid Heydarnia
% Date: June 12, 2025
% Last Modified: June 12, 2025


%   This code is an interpretation and adaptation of concepts and calculations
%   from the 'forces and torques' block from the original MegAWES repository 
%   (Copyright 2021 Delft University of Technology), which is also licensed under
%   the Apache License, Version 2.0.

function [F_a_B, M_a_B, AeroAngels] = aero_MegAWES(Omega_OB_B, v_a_W, alpha, beta, rho_factor, deflection_act, params_sym)
    % AERO_MEGAWES Calculates aerodynamic forces and moments for a multi-element kite.
    %
    % This function computes the total aerodynamic forces and moments acting on
    % a kite, modeled as a main wing, an elevator, and a rudder. It aggregates
    % the contributions from each component.
    %
    % Inputs:
    %   Omega_OB_B      - Vector (3x1), Angular velocity of the body frame 'B'
    %                     relative to the orientation frame 'O', expressed in body coordinates (rad/s).
    %   v_a_W           - Vector (3x1), Apparent wind velocity vector in the World frame (m/s).
    %   alpha           - Scalar, Angle of attack (rad).
    %   beta            - Scalar, Sideslip angle (rad).
    %   rho_factor      - Scalar, Factor for air density variation (e.g., due to altitude).
    %   deflection_act  - Vector (3x1), Actual control surface deflections (rad):
    %                     deflection_act(1) - Aileron deflection.
    %                     deflection_act(2) - Elevator deflection.
    %                     deflection_act(3) - Rudder deflection.
    %   params_sym      - struct, Symbolic parameters containing kite geometry
    %                     and environmental wind properties.
    %
    % Outputs:
    %   F_a_B           - Vector (3x1), Total aerodynamic force vector on the kite,
    %                     expressed in the body frame (N). [Fx; Fy; Fz]
    %   M_a_B           - Vector (3x1), Total aerodynamic moment vector on the kite,
    %                     expressed in the body frame (Nm). [Mx; My; Mz]
    %   AeroAngels      - struct, Contains aerodynamic angles for output:
    %                     AeroAngels.ai   - Angle of attack (deg).
    %                     AeroAngels.el   - Corrected elevator angle of attack (deg).
    %                     AeroAngels.ru   - Corrected rudder sideslip angle (deg).
    %                     AeroAngels.defl - Corrected aileron deflection (deg).

    % Extract kite geometric parameters
    b = params_sym.kite.b;     % Wing span
    c = params_sym.kite.c;     % Mean aerodynamic chord
    b_el = params_sym.kite.b_el; % Elevator span
    c_el = params_sym.kite.c_el; % Elevator chord
    x_el = params_sym.kite.x_el; % Longitudinal position of elevator
    
    b_vs = params_sym.kite.b_vs; % Vertical stabilizer (rudder) span
    c_vs = params_sym.kite.c_vs; % Vertical stabilizer (rudder) chord
    x_vs = params_sym.kite.x_vs; % Longitudinal position of vertical stabilizer
    
    % Calculate the magnitude of the apparent wind velocity.
    v_a = norm(v_a_W);
    
    % Extract individual control deflections
    delta_a = deflection_act(1); % Aileron deflection
    delta_e = deflection_act(2); % Elevator deflection
    delta_r = deflection_act(3); % Rudder deflection
    
    % --- Aileron and Main Wing (including Fluid-Structure Interaction effects) ---
    [out_Wing, defl_corrected] = Wing_Aero(v_a, alpha, beta, Omega_OB_B, delta_a, rho_factor, params_sym);
    
    % --- Elevator Aerodynamics ---
    % Calculate damping contribution to effective angle of attack for elevator.
    % This accounts for the pitch rate (Omega_OB_B(1)) effect at the elevator's longitudinal position.
    Alpha_damp = atan2(((x_el + c_el * 0.25) * Omega_OB_B(1)), (v_a * cos(alpha) * cos(beta)));
    % Calculate the corrected angle of attack for the elevator in degrees.
    Alpha_corrected = (180 / pi) * (alpha + Alpha_damp + delta_e);
    
    % Polynomial terms for elevator aerodynamic coefficients (CL and CD)
    Poly_AoA = [Alpha_corrected^2; Alpha_corrected; 1];
    % Coefficients for elevator lift coefficient (CL_elevator)
    e_cL = [-4.5336e-05, 0.064178, 0.0019657];
    % Coefficients for elevator drag coefficient (CD_elevator)
    e_cD = [0.0016042, -0.00015091, 0.010216];
    
    % Calculate CL and CD for the elevator using polynomial interpolation
    CL_elevator = e_cL * Poly_AoA;
    CD_elevator = e_cD * Poly_AoA;
    
    % Prepare input vector for simpleAero_Elevator
    input_elevator = [v_a; alpha; beta; CL_elevator; CD_elevator];
    out_elevator = simpleAero_Elevator(input_elevator, rho_factor, params_sym);
    
    % --- Rudder Aerodynamics ---
    % Calculate damping contribution to effective sideslip angle for rudder.
    % This accounts for the yaw rate (Omega_OB_B(3)) effect at the rudder's longitudinal position.
    Beta_damp = -1 * atan2(((x_vs + c_vs * 0.25) * Omega_OB_B(3)), (v_a * cos(alpha) * cos(beta)));
    % Calculate the corrected sideslip angle for the rudder in degrees.
    % Note: Using x_el and c_el for rudder damping calculation might be a typo,
    % typically x_vs and c_vs (vertical stabilizer properties) would be used here.
    Beta_corrected = (180 / pi) * (beta + Beta_damp + delta_r);
    
    % Polynomial terms for rudder aerodynamic coefficients (CL and CD)
    Poly_SS = [Beta_corrected^2; Beta_corrected; 1];
    % Coefficients for rudder lift coefficient (CL_rudder)
    r_cL = [-2.7494e-05, 0.038884, 0.0011889];
    % Coefficients for rudder drag coefficient (CD_rudder)
    r_cD = [0.0039814, -0.00038564, 0.018132];
    
    % Calculate CL and CD for the rudder using polynomial interpolation
    CL_rudder = r_cL * Poly_SS;
    CD_rudder = r_cD * Poly_SS;
    
    % Prepare input vector for simpleAero_Rudder
    input_rudder = [v_a; alpha; beta; CL_rudder; CD_rudder];
    out_Rudder = simpleAero_Rudder(input_rudder, rho_factor, params_sym);
    
    % --- Total Aerodynamic Forces and Moments ---
    % Sum the forces and moments from the wing, elevator, and rudder.
    out_total = out_Wing + out_elevator + out_Rudder;
    
    % Assign the combined forces and moments to the output variables.
    F_a_B = out_total(1:3); % Total aerodynamic forces in body frame
    M_a_B = out_total(4:6); % Total aerodynamic moments in body frame

    % Populate the AeroAngels output structure with relevant angles.
    AeroAngels.ai   = alpha * 180 / pi; % Angle of attack in degrees
    AeroAngels.el   = Alpha_corrected;  % Corrected elevator angle of attack in degrees
    AeroAngels.ru   = Beta_corrected;   % Corrected rudder sideslip angle in degrees
    AeroAngels.defl = defl_corrected;   % Corrected aileron deflection in degrees
end

%---
% Sub-functions for aero_MegAWES
%---

function [out, defl_corrected] = Wing_Aero(v_a, alpha, beta, Omega_OB_B, Deflection_act_1, rho_factor, params_sym)
    % WING_AERO Calculates aerodynamic forces and moments for the main wing, including aileron effects.
    %
    % This function determines the forces and moments generated by the main wing
    % and ailerons, incorporating static aerodynamic coefficients and a damping
    % effect due to roll rate.
    %
    % Inputs:
    %   v_a             - Scalar, Magnitude of the apparent wind velocity (m/s).
    %   alpha           - Scalar, Angle of attack (rad).
    %   beta            - Scalar, Sideslip angle (rad). (Note: beta is used in damping calc, but not in static wing coefficients directly).
    %   Omega_OB_B      - Vector (3x1), Angular velocity of the body frame (rad/s).
    %   Deflection_act_1 - Scalar, Raw aileron deflection (rad).
    %   rho_factor      - Scalar, Factor for air density variation.
    %   params_sym      - struct, Symbolic parameters containing kite geometry.
    %
    % Outputs:
    %   out             - Vector (6x1), Forces and moments for the wing in body frame [Fx; Fy; Fz; Mx; My; Mz].
    %   defl_corrected  - Scalar, Corrected aileron deflection due to roll damping (rad).

    % Extract wing span
    b = params_sym.kite.b;
    
    % Convert angle of attack to degrees for polynomial evaluation
    alpha_deg = (180 / pi) * alpha;
    % Polynomial terms for angle of attack
    Poly_AoA = [alpha_deg^2; alpha_deg; 1];

    % Coefficients for static wing lift coefficient (CL_Wing_Static)
    ai_cL = [-0.00019837, 0.1033, 1.4277];
    % Coefficients for static wing drag coefficient (CD_Wing_Static)
    ai_cD = [0.0012436, 0.02186, 0.12672];
    % Coefficients for static wing pitching moment coefficient (CM_Wing_Static)
    ai_cM = [-8.4578e-06, -0.0093078, -0.36671];
    
    % Calculate static CL, CD, and CM for the wing using polynomial interpolation
    CL_Wing_Static = ai_cL * Poly_AoA;
    CD_Wing_Static = ai_cD * Poly_AoA;
    CM_Wing_Static = ai_cM * Poly_AoA;

    % Calculate damping contribution to aileron effectiveness due to roll rate.
    % The `0.5236` corresponds to 30 degrees (pi/6 radians), an estimation for maximum deflection.
    defl_damp = (atan2(((b * 0.75) * Omega_OB_B(1)), (v_a * cos(alpha) * cos(beta)))) / (0.5236);
    
    % Calculate the corrected aileron deflection.
    defl_corrected = Deflection_act_1 - defl_damp;
    
    % Calculate the static rolling moment coefficient due to the corrected aileron deflection.
    CR_Wing_Static = params_sym.kite.wing_cRoll_Static * defl_corrected;
    
    % Prepare input vector for simpleAero_Wing
    input_aileron = [v_a; alpha; CL_Wing_Static; CD_Wing_Static; CM_Wing_Static; CR_Wing_Static];
    
    % Call the core wing aerodynamic calculation function
    out = simpleAero_Wing(input_aileron, rho_factor, params_sym);
end

%---

function out_el = simpleAero_Elevator(input, rho_factor, params_sym)
    % SIMPLEAERO_ELEVATOR Calculates aerodynamic forces and moments for the elevator.
    %
    % This function computes the lift and drag forces on the elevator based on
    % its calculated lift and drag coefficients, then transforms these into
    % body-frame forces and moments.
    %
    % Inputs:
    %   input        - Vector (5x1), Contains:
    %                  input(1) - Apparent wind speed magnitude (v).
    %                  input(2) - Angle of attack (alpha, for force projection).
    %                  input(3) - Sideslip angle (beta, for force projection).
    %                  input(4) - Calculated lift coefficient (CL_elevator).
    %                  input(5) - Calculated drag coefficient (CD_elevator).
    %   rho_factor   - Scalar, Factor for air density variation.
    %   params_sym   - struct, Symbolic parameters containing kite geometry:
    %                  params_sym.env.wind.rho - Reference air density.
    %                  params_sym.kite.s_el    - Elevator reference area.
    %                  params_sym.kite.x_el    - Longitudinal position of elevator.
    %                  params_sym.kite.c_el    - Elevator chord.
    %
    % Outputs:
    %   out_el       - Vector (6x1), Forces and moments for the elevator in body frame [Fx; Fy; Fz; Mx; My; Mz].

    % Calculate actual air density
    rho = params_sym.env.wind.rho * rho_factor;
    
    % Extract inputs
    v     = input(1); % Apparent wind speed magnitude
    alpha = input(2); % Angle of attack
    beta  = input(3); % Sideslip angle
    
    % Calculate dynamic pressure
    q = 0.5 * rho * v^2;

    % Calculate lift and drag forces on the elevator (sign conventions need careful checking)
    % Negative sign for F_L_el and F_D_el likely means these are defined as positive in the
    % direction opposite to the conventional lift/drag in an aerodynamic frame,
    % or it relates to the body frame transformation later.
    F_L_el = -input(4) * q * params_sym.kite.s_el; % Elevator lift force
    F_D_el = -input(5) * q * params_sym.kite.s_el; % Elevator drag force

    % Transform forces from an aerodynamic frame (or component-aligned frame) to the body frame.
    % Assuming elevator lift acts primarily in Z-body and drag in X-body, then projected.
    F_x_s = F_D_el * cos(alpha) * cos(beta) - F_L_el * sin(alpha) * cos(beta);
    F_y_s = 0; % Assuming no side force from elevator
    F_z_s = F_L_el * cos(alpha); % This projection seems simplified; careful with sign.

    % Calculate moments generated by the elevator in the body frame.
    % Assuming moments are about the center of gravity (COG) or reference point.
    % Only pitching moment (My) is considered, due to lift force.
    M_x_s = 0; % No roll moment from elevator (ideally)
    M_y_s = F_L_el * (params_sym.kite.x_el + params_sym.kite.c_el * 0.25); % Pitching moment
    M_z_s = 0; % No yaw moment from elevator (ideally)

    % Combine forces and moments into a single output vector
    out_el = [F_x_s, F_y_s, F_z_s, M_x_s, M_y_s, M_z_s]';
end

%---

function out_ru = simpleAero_Rudder(input, rho_factor, params_sym)
    % SIMPLEAERO_RUDDER Calculates aerodynamic forces and moments for the rudder (vertical stabilizer).
    %
    % This function computes the lift and drag forces on the rudder based on
    % its calculated coefficients, then transforms these into body-frame forces
    % and moments.
    %
    % Inputs:
    %   input        - Vector (5x1), Contains:
    %                  input(1) - Apparent wind speed magnitude (v).
    %                  input(2) - Angle of attack (alpha, for force projection).
    %                  input(3) - Sideslip angle (beta, for force projection).
    %                  input(4) - Calculated lift coefficient (CL_rudder).
    %                  input(5) - Calculated drag coefficient (CD_rudder).
    %   rho_factor   - Scalar, Factor for air density variation.
    %   params_sym   - struct, Symbolic parameters containing kite geometry:
    %                  params_sym.env.wind.rho - Reference air density.
    %                  params_sym.kite.s_vs    - Vertical stabilizer reference area.
    %                  params_sym.kite.x_vs    - Longitudinal position of vertical stabilizer.
    %                  params_sym.kite.c_vs    - Vertical stabilizer chord.
    %                  params_sym.kite.b_vs    - Vertical stabilizer span.
    %
    % Outputs:
    %   out_ru       - Vector (6x1), Forces and moments for the rudder in body frame [Fx; Fy; Fz; Mx; My; Mz].

    % Calculate actual air density
    rho = params_sym.env.wind.rho * rho_factor;
    
    % Extract inputs
    v     = input(1); % Apparent wind speed magnitude
    alpha = input(2); % Angle of attack
    beta  = input(3); % Sideslip angle

    % Calculate dynamic pressure
    q = 0.5 * rho * v^2;
    
    % Calculate lift and drag forces on the rudder (sign conventions need careful checking).
    % The factor of 2 for F_L_vs suggests modeling of a vertical stabilizer, potentially two surfaces
    % or specific coefficient definition.
    F_L_vs = 2 * input(4) * q * params_sym.kite.s_vs; % Rudder lift force (side force)
    F_D_vs = -2 * input(5) * q * params_sym.kite.s_vs; % Rudder drag force

    % Transform forces from an aerodynamic frame (or component-aligned frame) to the body frame.
    % Assuming rudder lift acts primarily in Y-body (side force) and drag in X-body, then projected.
    F_x_s = F_D_vs * cos(alpha) * cos(beta); % Drag contribution to X-force
    F_y_s = -F_L_vs * cos(beta); % Lift contribution to Y-force (side force)
    F_z_s = 0; % No Z-force from rudder (ideally)

    % Calculate moments generated by the rudder in the body frame.
    % Assuming moments are about the center of gravity (COG) or reference point.
    M_x_s = F_L_vs * (params_sym.kite.b_vs * 0.5); % Roll moment due to side force
    M_y_s = 0; % No pitch moment from rudder (ideally)
    M_z_s = F_L_vs * (params_sym.kite.x_vs + params_sym.kite.c_vs * 0.25); % Yaw moment due to side force

    % Combine forces and moments into a single output vector
    out_ru = [F_x_s, F_y_s, F_z_s, M_x_s, M_y_s, M_z_s]';
end

%---

function out_w = simpleAero_Wing(input, rho_factor, params_sym)
    % SIMPLEAERO_WING Calculates aerodynamic forces and moments for the main wing.
    %
    % This function computes the lift and drag forces, and rolling and pitching moments
    % on the main wing, then transforms these into body-frame forces and moments.
    % It uses input coefficients (CL, CD, CM_pitch, CL_roll) and dynamic pressure.
    %
    % Inputs:
    %   input        - Vector (6x1), Contains:
    %                  input(1) - Apparent wind speed magnitude (v).
    %                  input(2) - Angle of attack (alpha).
    %                  input(3) - Calculated lift coefficient (CL_Wing_Static).
    %                  input(4) - Calculated drag coefficient (CD_Wing_Static).
    %                  input(5) - Calculated pitching moment coefficient (CM_Wing_Static).
    %                  input(6) - Calculated rolling moment coefficient (CR_Wing_Static).
    %   rho_factor   - Scalar, Factor for air density variation.
    %   params_sym   - struct, Symbolic parameters containing kite geometry:
    %                  params_sym.kite.b       - Wing span.
    %                  params_sym.kite.c       - Mean aerodynamic chord.
    %                  params_sym.kite.s_ref   - Reference wing area.
    %                  params_sym.env.wind.rho - Reference air density.
    %                  params_sym.kite.cog(1)  - X-coordinate of center of gravity (for moment arm).
    %
    % Outputs:
    %   out_w        - Vector (6x1), Forces and moments for the wing in body frame [Fx; Fy; Fz; Mx; My; Mz].

    % Extract inputs
    v     = input(1); % Apparent wind speed magnitude
    alpha = input(2); % Angle of attack

    % Extract kite geometry and air density
    b     = params_sym.kite.b;     % Wing span
    c     = params_sym.kite.c;     % Mean aerodynamic chord
    s_ref = params_sym.kite.s_ref; % Reference wing area
    rho   = params_sym.env.wind.rho * rho_factor; % Actual air density
    
    % Calculate dynamic pressure
    q = 0.5 * rho * v^2;

    % Calculate main lift and drag forces
    F_L_main = -input(3) * q * s_ref; % Main wing lift force (negative sign likely due to body frame convention)
    F_D_main = input(4) * q * s_ref;  % Main wing drag force

    % Calculate pitching moment (My) and rolling moment (Mx)
    M_y_mw = input(5) * q * s_ref * c + ... % Static pitching moment from CM_Wing_Static
             (0.25 * c + params_sym.kite.cog(1)) * F_L_main; % Pitching moment due to lift at COG
             
    c_roll_mw = input(6); % Rolling moment coefficient
    M_x_mw = (c_roll_mw) * q * s_ref * b; % Rolling moment

    % Transform lift and drag forces into body frame X and Z components.
    % Assumes lift acts along -Z_body and drag along X_body in the aerodynamic frame.
    F_x_mw = -F_D_main * cos(alpha) - F_L_main * sin(alpha);
    F_z_mw = F_L_main * cos(alpha);
    
    % Combine forces and moments into a single output vector.
    % Note: Fy (side force) and Mz (yaw moment) are assumed zero for the main wing in this model.
    out_w = [F_x_mw, 0, F_z_mw, M_x_mw, M_y_mw, 0]';
end