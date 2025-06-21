function [CX, CY, CZ, Cl, Cm, Cn] = set_borne_aero_params_num(aero_model)
% SET_BORNE_AERO_PARAMS_NUM Sets numerical aerodynamic coefficients for various models.
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
%   This function fills in numerical structures with aerodynamic coefficients
%   (CX, CY, CZ, Cl, Cm, Cn) based on the specified `aero_model`. These coefficients
%   are typically used in simplified aerodynamic models.
%
%   The 'ALM', 'VLM', and 'CFD' cases refer to different sets of coefficients,
%   likely derived from:
%   - 'ALM': Actuator Line method.
%            More info: Crismer, J.-B., et al. Large-Eddy Simulation of airborne wind energy systems wakes.
%                        in Journal of Physics: Conference Series. 2023. IOP Publishing. DOI 10.1088/1742-6596/2505/1/012036  

%   - 'VLM': Vortex Lattice Method
%            More info: Drela, M.a.H.Y. {AVL} 3.40b. 2024; Available from: http://web.mit.edu/drela/Public/web/avl/.
%
%   - 'CFD': Computational Fluid Dynamics.
%            More info: Pynaert, N., et al., Wing Deformation of an Airborne Wind Energy System in Crosswind Flight 
%                       Using High-Fidelity Fluid–Structure Interaction. Energies, 2023. 16(2): p. 602

% Inputs:
%   aero_model - string, Specifies which set of aerodynamic coefficients to load.
%                Valid options are: 'ALM', 'VLM', 'CFD'.
%
% Outputs:
%   CX         - struct, Coefficients for the X-body axis force coefficient. 
%   CY         - struct, Coefficients for the Y-body axis force coefficient.
%   CZ         - struct, Coefficients for the Z-body axis force coefficient.
%   Cl         - struct, Coefficients for the Rolling moment coefficient.
%   Cm         - struct, Coefficients for the Pitching moment coefficient.
%   Cn         - struct, Coefficients for the Yawing moment coefficient.
%

%
% See also: system_vars_params

% Use a switch-case statement to select the appropriate set of
% aerodynamic coefficients based on the input `aero_model` string.
    switch aero_model

        case 'ALM'
            % Aerodynamic coefficients derived from an Aerodynamic line method (ALM).

            % CX (Drag Coefficient) Components
            CX = struct();
            CX.zero = [-0.044491];                   % Baseline drag coefficient.
            CX.alpha = [0.643491, 4.312588];         % Coefficients for angle of attack (alpha) dependency.
            CX.q = [-0.435813, 3.576224, 9.319524];   % Coefficients for pitch rate (q) dependency.
            CX.deltae = [-0.032819, 0.382206, 0.336074]; % Coefficients for elevator deflection (delta_e) dependency.

            % CY (Side Force Coefficient) Components
            CY = struct();
            CY.beta = [-0.176968, -0.005307, 0.07225];     % Coefficients for sideslip angle (beta) dependency.
            CY.p = [-0.004113, -0.063994, -0.191728];     % Coefficients for roll rate (p) dependency.
            CY.r = [0.075628, -0.001207, -0.002044];     % Coefficients for yaw rate (r) dependency.
            CY.deltaa = [-0.000378, -0.000256, 0.002923]; % Coefficients for aileron deflection (delta_a) dependency.
            CY.deltar = [0.175971, 0.0086, -0.147732];     % Coefficients for rudder deflection (delta_r) dependency.

            % CZ (Lift Coefficient) Components
            CZ = struct();
            CZ.zero = [-0.975917];                   % Baseline lift coefficient.
            CZ.alpha = [-4.482273, 2.927043];         % Coefficients for angle of attack (alpha) dependency.
            CZ.q = [-2.862521, 4.477878, 49.332075];  % Coefficients for pitch rate (q) dependency.
            CZ.deltae = [-0.47211, 0.243975, 1.718238]; % Coefficients for elevator deflection (delta_e) dependency.

            % Cl (Rolling Moment Coefficient) Components
            Cl = struct();
            Cl.beta = [-0.007379, -6.2e-05, 0.00364];     % Coefficients for sideslip angle (beta) dependency.
            Cl.p = [-0.587414, 1.808228, 8.431355];     % Coefficients for roll rate (p) dependency.
            Cl.r = [0.244036, 0.613158, -0.716677];     % Coefficients for yaw rate (r) dependency.
            Cl.deltaa = [-0.370799, 0.001332, 0.177379]; % Coefficients for aileron deflection (delta_a) dependency.
            Cl.deltar = [0.007332, 0.000383, -0.006009]; % Coefficients for rudder deflection (delta_r) dependency.

            % Cm (Pitching Moment Coefficient) Components
            Cm = struct();
            Cm.zero = [0.057198];                    % Baseline pitching moment coefficient.
            Cm.alpha = [-0.139756, -0.376294];        % Coefficients for angle of attack (alpha) dependency.
            Cm.q = [-6.009464, -2.451568, -26.813118]; % Coefficients for pitch rate (q) dependency.
            Cm.deltae = [-1.201313, -0.016661, 0.043642]; % Coefficients for elevator deflection (delta_e) dependency.

            % Cn (Yawing Moment Coefficient) Components
            Cn = struct();
            Cn.beta = [0.040269, 0.00141, -0.012721];     % Coefficients for sideslip angle (beta) dependency.
            Cn.p = [-0.067101, -0.999124, -1.066629];     % Coefficients for roll rate (p) dependency.
            Cn.r = [-0.032551, 0.077463, 0.330239];     % Coefficients for yaw rate (r) dependency.
            Cn.deltaa = [0.004239, -0.338661, -0.009433]; % Coefficients for aileron deflection (delta_a) dependency.
            Cn.deltar = [-0.040178, -0.00205, 0.031991]; % Coefficients for rudder deflection (delta_r) dependency.

        case 'VLM'
            % Aerodynamic coefficients derived from a Vortex Lattice Method (VLM).

            % CX Components
            CX = struct();
            CX.zero = [-0.046];
            CX.alpha = [0.5329, 3.6178];
            CX.q = [-0.1689, 3.1142, -0.3229];
            CX.deltae = [-0.0203, 0.2281, 0.0541];

            % CY Components
            CY = struct();
            CY.beta = [-0.2056, -0.1529, -0.3609];
            CY.p = [0.0588, 0.3069, -0.0109];
            CY.r = [0.0869, 0.0271, -0.0541];
            CY.deltaa = [0.0064, -0.0365, -0.0022];
            CY.deltar = [0.1801, 0.0196, -0.1724];

            % CZ Components
            CZ = struct();
            CZ.zero = [-0.8781];
            CZ.alpha = [-4.7042, 0.0335];
            CZ.q = [-5.9365, -0.7263, 2.4422];
            CZ.deltae = [-0.4867, -0.007, 0.4642];

            % Cl (Rolling Moment Coefficient) Components
            Cl = struct();
            Cl.beta = [-0.0101, -0.1834, 0.0023];
            Cl.p = [-0.4888, -0.027, 0.092];
            Cl.r = [0.1966, 0.5629, -0.0498];
            Cl.deltaa = [-0.1972, 0.0574, 0.1674];
            Cl.deltar = [0.0077, -0.0091, -0.0092];

            % Cm (Pitching Moment Coefficient) Components
            Cm = struct();
            Cm.zero = [-0.065];
            Cm.alpha = [-0.3306, 0.2245];
            Cm.q = [-7.7531, -0.003, 3.8925];
            Cm.deltae = [-1.1885, -0.0007, 1.1612];

            % Cn (Yawing Moment Coefficient) Components
            Cn = struct();
            Cn.beta = [0.0385, 0.0001, -0.0441];
            Cn.p = [-0.0597, -0.7602, 0.0691];
            Cn.r = [-0.0372, -0.0291, -0.2164];
            Cn.deltaa = [0.0054, -0.0425, 0.0354];
            Cn.deltar = [-0.0404, -0.0031, 0.0385];

        case 'CFD'
            % Aerodynamic coefficients derived from Computational Fluid Dynamics (CFD).
            % Note: Some coefficients are commented out, indicating they may
            % not be used or available for this specific CFD dataset.

            % CX Components
            CX = struct();
            CX.zero = [-0.1164];
            CX.alpha = [0.4564, 2.3044];
            % CX.beta = [0.0279, 0.0414, 0.8307]; 
            % CX.p = [0.0342, 0.1529, -1.8588];   
            CX.q = [-0.4645, 8.5417, -10.8181];
            % CX.r = [-0.0006, 0.0519, 0.4025];   
            % CX.deltaa = [-0.0168, 0.0733, 1.3335]; 
            CX.deltae = [0.0002, -0.0182, 0.41];
            % CX.deltar = [-0.0173, -0.015, -0.2922]; 

            % CY Components
            CY = struct();
            % CY.'0' = [-0.0];                  
            % CY.alpha = [0.0002, 0.0013];      
            CY.beta = [-0.274, 0.1664, 0.8803];
            CY.p = [0.0198, -0.2312, -0.315];
            % CY.q = [0.0007, -0.001, 0.0799];   
            CY.r = [0.0911, -0.0267, -0.4982];
            CY.deltaa = [0.0063, 0.0119, -0.0754];
            % CY.deltae = [0.0001, -0.0012, -0.0216]; 
            CY.deltar = [0.2259, -0.1198, 0.1955];

            % CZ Components
            CZ = struct();
            CZ.zero = [-0.9245];
            CZ.alpha = [-3.7205, 4.7972];
            % CZ.beta = [0.1123, -0.125, -5.0971]; 
            % CZ.p = [0.1387, 0.1685, -27.9934]; 
            CZ.q = [-5.6405, 60.997, 240.6406];
            % CZ.r = [0.0067, 0.1349, -4.4412];   
            % CZ.deltaa = [0.0638, -1.8662, -26.6776]; 
            CZ.deltae = [-0.4897, 0.2366, 3.4195];
            % CZ.deltar = [0.0044, 0.0123, -0.2717]; 

            % Cl (Rolling Moment Coefficient) Components
            Cl = struct();
            % Cl.'0' = [0.0];                   
            % Cl.alpha = [0.0002, 0.0001];      
            Cl.beta = [0.0344, -0.1786, -2.6711];
            Cl.p = [-0.4052, 0.4109, -0.5721];
            % Cl.q = [0.018, 0.0258, -2.1828];   
            Cl.r = [0.1802, 0.5792, -0.0129];
            Cl.deltaa = [-0.0941, -0.1921, -0.2034];
            % Cl.deltae = [0.0, -0.0063, -0.0912]; 
            Cl.deltar = [0.0106, -0.0214, -0.0874];

            % Cm (Pitching Moment Coefficient) Components
            Cm = struct();
            Cm.zero = [0.0279];
            Cm.alpha = [-0.5307, -0.9786];
            % Cm.beta = [-0.0184, 0.7392, 8.2241]; 
            % Cm.p = [0.0008, -0.1007, -0.0845];   
            Cm.q = [-8.0446, 1.1837, -20.8571];
            % Cm.r = [-0.0021, -0.2081, -2.4176]; 
            % Cm.deltaa = [0.0177, 0.9504, 4.4178]; 
            Cm.deltae = [-1.2524, -0.092, 11.6916];
            % Cm.deltar = [0.0165, 0.0416, 0.0795]; 

            % Cn (Yawing Moment Coefficient) Components
            Cn = struct();
            % Cn.'0' = [-0.0];                   
            % Cn.alpha = [-0.0, 0.0004];         
            Cn.beta = [0.0682, 0.0048, -0.1193];
            Cn.p = [-0.0412, -0.4284, -1.0241];
            % Cn.q = [-0.0007, 0.0072, 0.0489];   
            Cn.r = [-0.0555, 0.0316, 0.1057];
            Cn.deltaa = [0.0234, -0.0113, -0.6566];
            % Cn.deltae = [-0.0, -0.0001, 0.0014]; 
            Cn.deltar = [-0.0509, 0.0287, -0.0572];
    end % End of switch aero_model

end % End of function set_borne_aero_params_num