function [F_a_B, M_a_B, AeroAngels] = aero(Omega_OB_B, V_a_W, alpha, beta, Deflection_act, F_T_total, Params)
    
    
    b = Params.DE2019.b;
    c = Params.DE2019.c;
    V_a_O = transformFromWtoO(Params.env.windDirection_rad,V_a_W);
    V_a = norm(V_a_W);
    inputs = struct();
    inputs.zero = 1;
    inputs.alpha = alpha;
    inputs.beta = beta;
    inputs.p = 0.5 * (1/V_a) * b * Omega_OB_B(1);
    inputs.q = 0.5 * (1/V_a) * c * Omega_OB_B(2);
    inputs.r = 0.5 * (1/V_a) * b * Omega_OB_B(3);
    inputs.deltaa = Deflection_act(1);
    inputs.deltae = Deflection_act(2);
    inputs.deltar = Deflection_act(3);
    
    switch Params.sim.aero_model
    
        case 'ALM'
             % CX
            stab_derivs.CX = struct();
            stab_derivs.CX.zero = [-0.044491];
            stab_derivs.CX.alpha = [0.643491, 4.312588];
            stab_derivs.CX.q = [-0.435813, 3.576224, 9.319524];
            stab_derivs.CX.deltae = [-0.032819, 0.382206, 0.336074];
            
            % CY
            stab_derivs.CY = struct();
            stab_derivs.CY.beta = [-0.176968, -0.005307, 0.07225];
            stab_derivs.CY.p = [-0.004113, -0.063994, -0.191728];
            stab_derivs.CY.r = [0.075628, -0.001207, -0.002044];
            stab_derivs.CY.deltaa = [-0.000378, -0.000256, 0.002923];
            stab_derivs.CY.deltar = [0.175971, 0.0086, -0.147732];
            
            % CZ
            stab_derivs.CZ = struct();
            stab_derivs.CZ.zero = [-0.975917];
            stab_derivs.CZ.alpha = [-4.482273, 2.927043];
            stab_derivs.CZ.q = [-2.862521, 4.477878, 49.332075];
            stab_derivs.CZ.deltae = [-0.47211, 0.243975, 1.718238];
            
            % Cl
            stab_derivs.Cl = struct();
            stab_derivs.Cl.beta = [-0.007379, -6.2e-05, 0.00364];
            stab_derivs.Cl.p = [-0.587414, 1.808228, 8.431355];
            stab_derivs.Cl.r = [0.244036, 0.613158, -0.716677];
            stab_derivs.Cl.deltaa = [-0.370799, 0.001332, 0.177379];
            stab_derivs.Cl.deltar = [0.007332, 0.000383, -0.006009];
            
            % Cm
            stab_derivs.Cm = struct();
            stab_derivs.Cm.zero = [0.057198];
            stab_derivs.Cm.alpha = [-0.139756, -0.376294];
            stab_derivs.Cm.q = [-6.009464, -2.451568, -26.813118];
            stab_derivs.Cm.deltae = [-1.201313, -0.016661, 0.043642];
            
            % Cn
            stab_derivs.Cn = struct();
            stab_derivs.Cn.beta = [0.040269, 0.00141, -0.012721];
            stab_derivs.Cn.p = [-0.067101, -0.999124, -1.066629];
            stab_derivs.Cn.r = [-0.032551, 0.077463, 0.330239];
            stab_derivs.Cn.deltaa = [0.004239, -0.338661, -0.009433];
            stab_derivs.Cn.deltar = [-0.040178, -0.00205, 0.031991];
            [F_a_B, M_a_B, AeroAngels] = BORNE_Aerodynamics(V_a_O, inputs, stab_derivs, F_T_total, Params);
    
        case 'VLM'
            % CX
            stab_derivs.CX = struct();
            stab_derivs.CX.zero = [-0.046];
            stab_derivs.CX.alpha = [0.5329, 3.6178];
            stab_derivs.CX.q = [-0.1689, 3.1142, -0.3229];
            stab_derivs.CX.deltae = [-0.0203, 0.2281, 0.0541];
            
            % CY
            stab_derivs.CY = struct();
            stab_derivs.CY.beta = [-0.2056, -0.1529, -0.3609];
            stab_derivs.CY.p = [0.0588, 0.3069, -0.0109];
            stab_derivs.CY.r = [0.0869, 0.0271, -0.0541];
            stab_derivs.CY.deltaa = [0.0064, -0.0365, -0.0022];
            stab_derivs.CY.deltar = [0.1801, 0.0196, -0.1724];
            
            % CZ
            stab_derivs.CZ = struct();
            stab_derivs.CZ.zero = [-0.8781];
            stab_derivs.CZ.alpha = [-4.7042, 0.0335];
            stab_derivs.CZ.q = [-5.9365, -0.7263, 2.4422];
            stab_derivs.CZ.deltae = [-0.4867, -0.007, 0.4642];
            
            % Cl
            stab_derivs.Cl = struct();
            stab_derivs.Cl.beta = [-0.0101, -0.1834, 0.0023];
            stab_derivs.Cl.p = [-0.4888, -0.027, 0.092];
            stab_derivs.Cl.r = [0.1966, 0.5629, -0.0498];
            stab_derivs.Cl.deltaa = [-0.1972, 0.0574, 0.1674];
            stab_derivs.Cl.deltar = [0.0077, -0.0091, -0.0092];
            
            % Cm
            stab_derivs.Cm = struct();
            stab_derivs.Cm.zero = [-0.065];
            stab_derivs.Cm.alpha = [-0.3306, 0.2245];
            stab_derivs.Cm.q = [-7.7531, -0.003, 3.8925];
            stab_derivs.Cm.deltae = [-1.1885, -0.0007, 1.1612];
            
            % Cn
            stab_derivs.Cn = struct();
            stab_derivs.Cn.beta = [0.0385, 0.0001, -0.0441];
            stab_derivs.Cn.p = [-0.0597, -0.7602, 0.0691];
            stab_derivs.Cn.r = [-0.0372, -0.0291, -0.2164];
            stab_derivs.Cn.deltaa = [0.0054, -0.0425, 0.0354];
            stab_derivs.Cn.deltar = [-0.0404, -0.0031, 0.0385];
            [F_a_B, M_a_B, AeroAngels] = BORNE_Aerodynamics(V_a_O, inputs, stab_derivs, F_T_total, Params);
    
        case 'CFD'
            % CX derivatives
            stab_derivs.CX = struct();
            stab_derivs.CX.zero = [-0.1164];
            stab_derivs.CX.alpha = [0.4564, 2.3044];
            % stab_derivs.CX.beta = [0.0279, 0.0414, 0.8307];
            % stab_derivs.CX.p = [0.0342, 0.1529, -1.8588];
            stab_derivs.CX.q = [-0.4645, 8.5417, -10.8181];
            % stab_derivs.CX.r = [-0.0006, 0.0519, 0.4025];
            % stab_derivs.CX.deltaa = [-0.0168, 0.0733, 1.3335];
            stab_derivs.CX.deltae = [0.0002, -0.0182, 0.41];
            % stab_derivs.CX.deltar = [-0.0173, -0.015, -0.2922];
            
            % CY derivatives
            stab_derivs.CY = struct();
            % stab_derivs.CY.'0' = [-0.0];
            % stab_derivs.CY.alpha = [0.0002, 0.0013];
            stab_derivs.CY.beta = [-0.274, 0.1664, 0.8803];
            stab_derivs.CY.p = [0.0198, -0.2312, -0.315];
            % stab_derivs.CY.q = [0.0007, -0.001, 0.0799];
            stab_derivs.CY.r = [0.0911, -0.0267, -0.4982];
            stab_derivs.CY.deltaa = [0.0063, 0.0119, -0.0754];
            % stab_derivs.CY.deltae = [0.0001, -0.0012, -0.0216];
            stab_derivs.CY.deltar = [0.2259, -0.1198, 0.1955];
            
            % CZ derivatives
            stab_derivs.CZ = struct();
            stab_derivs.CZ.zero = [-0.9245];
            stab_derivs.CZ.alpha = [-3.7205, 4.7972];
            % stab_derivs.CZ.beta = [0.1123, -0.125, -5.0971];
            % stab_derivs.CZ.p = [0.1387, 0.1685, -27.9934];
            stab_derivs.CZ.q = [-5.6405, 60.997, 240.6406];
            % stab_derivs.CZ.r = [0.0067, 0.1349, -4.4412];
            % stab_derivs.CZ.deltaa = [0.0638, -1.8662, -26.6776];
            stab_derivs.CZ.deltae = [-0.4897, 0.2366, 3.4195];
            % stab_derivs.CZ.deltar = [0.0044, 0.0123, -0.2717];
            
            % Cl derivatives
            stab_derivs.Cl = struct();
            % stab_derivs.Cl.'0' = [0.0];
            % stab_derivs.Cl.alpha = [0.0002, 0.0001];
            stab_derivs.Cl.beta = [0.0344, -0.1786, -2.6711];
            stab_derivs.Cl.p = [-0.4052, 0.4109, -0.5721];
            % stab_derivs.Cl.q = [0.018, 0.0258, -2.1828];
            stab_derivs.Cl.r = [0.1802, 0.5792, -0.0129];
            stab_derivs.Cl.deltaa = [-0.0941, -0.1921, -0.2034];
            % stab_derivs.Cl.deltae = [0.0, -0.0063, -0.0912];
            stab_derivs.Cl.deltar = [0.0106, -0.0214, -0.0874];
            
            % Cm derivatives
            stab_derivs.Cm = struct();
            stab_derivs.Cm.zero = [0.0279];
            stab_derivs.Cm.alpha = [-0.5307, -0.9786];
            % stab_derivs.Cm.beta = [-0.0184, 0.7392, 8.2241];
            % stab_derivs.Cm.p = [0.0008, -0.1007, -0.0845];
            stab_derivs.Cm.q = [-8.0446, 1.1837, -20.8571];
            % stab_derivs.Cm.r = [-0.0021, -0.2081, -2.4176];
            % stab_derivs.Cm.deltaa = [0.0177, 0.9504, 4.4178];
            stab_derivs.Cm.deltae = [-1.2524, -0.092, 11.6916];
            % stab_derivs.Cm.deltar = [0.0165, 0.0416, 0.0795];
            
            % Cn derivatives
            stab_derivs.Cn = struct();
            % stab_derivs.Cn.'0' = [-0.0];
            % stab_derivs.Cn.alpha = [-0.0, 0.0004];
            stab_derivs.Cn.beta = [0.0682, 0.0048, -0.1193];
            stab_derivs.Cn.p = [-0.0412, -0.4284, -1.0241];
            % stab_derivs.Cn.q = [-0.0007, 0.0072, 0.0489];
            stab_derivs.Cn.r = [-0.0555, 0.0316, 0.1057];
            stab_derivs.Cn.deltaa = [0.0234, -0.0113, -0.6566];
            % stab_derivs.Cn.deltae = [-0.0, -0.0001, 0.0014];
            stab_derivs.Cn.deltar = [-0.0509, 0.0287, -0.0572];

            [F_a_B, M_a_B, AeroAngels] = BORNE_Aerodynamics(V_a_W, inputs, stab_derivs, F_T_total, Params);
        otherwise % xfoil -> MegAWES built in aerodynamic model
            [F_a_B, M_a_B, AeroAngels] = Aerodynamics(Omega_OB_B, V_a, alpha, beta, Deflection_act, F_T_total, Params);
    end

end
function [F_a_B, M_a_B, AeroAngels] = BORNE_Aerodynamics(V_a_O, inputs, stab_derivs, F_T_total, Params)
    S_ref = Params.DE2019.S_ref;
    b = Params.DE2019.b;
    c = Params.DE2019.c;
    aero_coeff = struct();
    fields_l0 = fieldnames(stab_derivs);
    for i = 1:length(fields_l0)
    field_0 = fields_l0{i}; % CX CY CZ Cl Cm Cn
    fields_l1 = fieldnames(stab_derivs.(field_0));
    aero_coeff.(field_0) = 0;
        for j = 1:length(fields_l1)
               field_1 = fields_l1{j}; % zero, alpha, beta, ...
               Alpha = []; 
               for len=1:length(stab_derivs.(field_0).(field_1))
                   Alpha = [Alpha; inputs.alpha^(len-1)];
               end
               aero_coeff.(field_0) = aero_coeff.(field_0) + (stab_derivs.(field_0).(field_1) * Alpha) * inputs.(field_1); 
        end
    end
%     q_eff = 0.5 * Params.env.rhos * (V_a_O' * V_a_O) * S_ref;
    Va = norm(V_a_O);
    q_eff = 0.5 * Params.env.rhos * Va^2 * S_ref;
    F_a_B = q_eff * [aero_coeff.CX; aero_coeff.CY; aero_coeff.CZ]; 
    M_a_B = q_eff * [b * aero_coeff.Cl; c * aero_coeff.Cm; b * aero_coeff.Cn];

     if Params.sim.ablation_study == 1
        F_a_B = F_a_B + F_T_total;
     end
    AeroAngels = struct();
    AeroAngels.ai = inputs.alpha * 180/pi;
    AeroAngels.ru = inputs.beta * 180/pi;
end



