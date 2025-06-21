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
% Email: omid.heydarnia@ugent.be
% Date: June 12, 2025
% Last Modified: June 12, 2025

function [F_aero_B, M_aero_B, AeroAngels] = aero_borne(omega_OB_B, v_a_W, alpha, beta, rho_factor, deflection_act, params_sym)
    % AERO_BORNE Calculates aerodynamic forces and moments on the airborne system.
    %
    % This function serves as a wrapper to prepare inputs for the core
    % BORNE_Aerodynamics function. It transforms the apparent wind velocity
    % from the World frame to a local aerodynamic frame and organizes all
    % necessary parameters and state variables.
    %
    % Inputs:
    %   omega_OB_B      - Vector (3x1), Angular velocity of the body frame 'B'
    %                     relative to the orientation frame 'O', expressed in body coordinates (rad/s).
    %                     [p; q; r] in aerospace convention.
    %   v_a_W           - Vector (3x1), Apparent wind velocity vector in the wind frame (m/s).
    %   alpha           - Scalar, Angle of attack (rad).
    %   beta            - Scalar, Sideslip angle (rad).
    %   rho_factor      - Scalar, Factor for air density variation (e.g., due to altitude).
    %   deflection_act  - Vector (3x1), Actual control surface deflections (rad):
    %                     deflection_act(1) - Aileron deflection.
    %                     deflection_act(2) - Elevator deflection.
    %                     deflection_act(3) - Rudder deflection.
    %   params_sym      - struct, Symbolic parameters.
    %
    % Outputs:
    %   F_aero_B        - Vector (3x1), Total aerodynamic force vector on the kite,
    %                     expressed in the body frame (N).
    %   M_aero_B        - Vector (3x1), Total aerodynamic moment vector on the kite,
    %                     expressed in the body frame (Nm).
    %   AeroAngels      - struct, Contains aerodynamic angles in degrees for output:
    %                     AeroAngels.ai   - Angle of attack (deg).
    %                     AeroAngels.ru   - Sideslip angle (deg).
    %                     AeroAngels.defl - Control deflections (deg).

    % Extract kite geometry parameters
    b = params_sym.kite.b; % Wing span
    c = params_sym.kite.c; % Mean aerodynamic chord


    % Calculate the magnitude of the apparent wind velocity.
    v_a = norm(v_a_W);

    % Assemble all necessary inputs for the aerodynamic model into a structure.
    % These inputs include aerodynamic angles, body rates (normalized), and control deflections.
    inputs = struct();
    inputs.zero = 1; % Placeholder/constant term for stability derivative summation
    inputs.alpha = alpha; % Angle of attack
    inputs.beta = beta;   % Sideslip angle
    % Normalized body rates (non-dimensional):
    inputs.p = 0.5 * (1/v_a) * b * omega_OB_B(1); % Normalized roll rate (p_bar)
    inputs.q = 0.5 * (1/v_a) * c * omega_OB_B(2); % Normalized pitch rate (q_bar)
    inputs.r = 0.5 * (1/v_a) * b * omega_OB_B(3); % Normalized yaw rate (r_bar)
    % Control surface deflections:
    inputs.delta_a = deflection_act(1); % Aileron deflection
    inputs.delta_e = deflection_act(2); % Elevator deflection
    inputs.delta_r = deflection_act(3); % Rudder deflection

    % Extract aerodynamic stability derivatives from parameters.
    % These structs contain coefficients for forces (CX, CY, CZ) and moments (Cl, Cm, Cn).
    stab_derivs = struct();
    stab_derivs.CX = params_sym.kite.CX;
    stab_derivs.CY = params_sym.kite.CY;
    stab_derivs.CZ = params_sym.kite.CZ;
    stab_derivs.Cl = params_sym.kite.Cl;
    stab_derivs.Cm = params_sym.kite.Cm;
    stab_derivs.Cn = params_sym.kite.Cn;

    % Call the core aerodynamic calculation function.
    [F_aero_B, M_aero_B, AeroAngels] = BORNE_Aerodynamics(v_a_W, inputs, stab_derivs, rho_factor, params_sym);
end


function [F_a_B, M_a_B, AeroAngels] = BORNE_Aerodynamics(v_a_W, inputs, stab_derivs, rho_factor, params_sym)
    % BORNE_AERODYNAMICS Calculates aerodynamic forces and moments based on stability derivatives.
    %
    % This function computes the aerodynamic coefficients for forces (CX, CY, CZ)
    % and moments (Cl, Cm, Cn) by summing contributions from various aerodynamic
    % inputs (angles, rates, deflections) based on a linear or polynomial
    % representation of stability derivatives. It then uses these coefficients
    % and dynamic pressure to calculate the total forces and moments.
    %
    % Inputs:
    %   v_a_W       - Vector (3x1), Apparent wind velocity vector in the wind frame (m/s).
    %   inputs      - struct, Contains all aerodynamic inputs:
    %                 inputs.zero, inputs.alpha, inputs.beta, inputs.p, inputs.q, inputs.r,
    %                 inputs.delta_a, inputs.delta_e, inputs.delta_r.
    %   stab_derivs - struct, Contains aerodynamic stability derivatives:
    %                 CX, CY, CZ, Cl, Cm, Cn (each of which is a struct).
    %   rho_factor  - Scalar, Factor for air density variation (e.g., due to altitude).
    %   params_sym  - struct, Symbolic parameters.
    %
    % Outputs:
    %   F_a_B       - Vector (3x1), Aerodynamic force vector in the body frame (N).
    %   M_a_B       - Vector (3x1), Aerodynamic moment vector in the body frame (Nm).
    %   AeroAngels  - struct, Contains aerodynamic angles in degrees for output:
    %                 AeroAngels.ai   - Angle of attack (deg).
    %                 AeroAngels.ru   - Sideslip angle (deg).
    %                 AeroAngels.defl - Control deflections (deg).

    % Extract kite geometry and environmental parameters
    s_ref = params_sym.kite.s_ref; % Reference wing area
    b     = params_sym.kite.b;     % Wing span
    c     = params_sym.kite.c;     % Mean aerodynamic chord
    rho   = params_sym.env.wind.rho * rho_factor; % Actual air density

    % Initialize structure to store calculated aerodynamic coefficients (CX, CY, CZ, Cl, Cm, Cn)
    aero_coeff = struct();

    % Iterate through each main aerodynamic coefficient (CX, CY, CZ, Cl, Cm, Cn)
    fields_l0 = fieldnames(stab_derivs); % e.g., {'CX', 'CY', 'CZ', 'Cl', 'Cm', 'Cn'}
    for i = 1:length(fields_l0)
        field_0 = fields_l0{i}; % Current coefficient name (e.g., 'CX')
        
        % Get the names of the input terms contributing to this coefficient (e.g., 'zero', 'alpha', 'beta', etc.)
        fields_l1 = fieldnames(stab_derivs.(field_0));
        aero_coeff.(field_0) = 0; % Initialize current coefficient to zero

        % Iterate through each contributing input term (e.g., 'zero', 'alpha', 'beta', 'p', 'q', 'r', 'delta_a', etc.)
        for j = 1:length(fields_l1)
            field_1 = fields_l1{j}; % Current input term name (e.g., 'alpha')
            Alpha = []; % Initialize matrix for alpha powers
            
            % Build a column vector of powers of alpha [1; alpha; alpha^2; ...]
            % This allows for polynomial dependence of derivatives on alpha.
            for len = 1:length(stab_derivs.(field_0).(field_1))
                Alpha = [Alpha; inputs.alpha^(len-1)];
            end
            
            % Calculate the contribution of the current input term to the total coefficient.
            % This is a sum of (derivative_coefficient * alpha_power) * input_value.
            % (stab_derivs.(field_0).(field_1).') is the row vector of coefficients for this input term.
            % Alpha is the column vector of alpha powers.
            % (stab_derivs.(field_0).(field_1).' * Alpha) performs the polynomial evaluation for the derivative.
            % The result is then multiplied by the actual input value (e.g., inputs.alpha, inputs.beta, etc.).
            aero_coeff.(field_0) = aero_coeff.(field_0) + (stab_derivs.(field_0).(field_1).' * Alpha) * inputs.(field_1); 
        end
    end

    % Calculate the effective dynamic pressure times reference area (q_eff).
    % This term converts non-dimensional aerodynamic coefficients into dimensional foces/moments.
    q_eff = 0.5 * rho * (v_a_W' * v_a_W) * s_ref;

    % Calculate total aerodynamic forces in the body frame.
    % CX, CY, CZ are body-axis force coefficients.
    F_a_B = q_eff * [aero_coeff.CX; aero_coeff.CY; aero_coeff.CZ]; 
    
    % Calculate total aerodynamic moments in the body frame.
    % Cl, Cm, Cn are body-axis moment coefficients, scaled by span (b) or chord (c).
    M_a_B = q_eff * [b * aero_coeff.Cl; c * aero_coeff.Cm; b * aero_coeff.Cn];

    % Populate AeroAngels struct with relevant angles converted to degrees for easier interpretation.
    AeroAngels = struct();
    AeroAngels.ai   = inputs.alpha * 180/pi; % Angle of attack in degrees
    AeroAngels.ru   = inputs.beta * 180/pi;   % Sideslip angle in degrees
    AeroAngels.defl = [inputs.delta_a; inputs.delta_e; inputs.delta_r] .* 180/pi; % Control deflections in degrees
end