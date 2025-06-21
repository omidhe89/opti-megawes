function save_flight_data(t, optimal_resutls, sol_stat, params_num)
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% See the accompanying LICENSE file for the full text of the license.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 10, 2025
% Last Modified: June 10, 2025

    
    baseWindSpeed = params_num.env.wind.base_velocity;
  
    
    if params_num.sim.path_type == 0
        if params_num.sim.ablation_study
            results = struct();
            results.optimal_sol = optimal_resutls;
            results.time = t;
            results.solstats = sol_stat;
            filename = strcat("./results/_ablation_0_", char(params_num.sim.aero_model), "_V_", num2str(baseWindSpeed), "Nw_", num2str(params_num.sim.winding_number) , ".mat");
            save(filename, "results");
        else
            results = struct();
            results.optimal_sol = optimal_resutls;
            results.time = t;
            results.solstats = sol_stat;
            filename = strcat("./results/_normal_0_", char(params_num.sim.aero_model), "_V_", num2str(baseWindSpeed), "_Nw_", num2str(params_num.sim.winding_number) , ".mat");
            save(filename, "results");
        end
    else
        if params_num.sim.ablation_study
            results = struct();
            results.optimal_sol = optimal_resutls;
            results.time = t;
            results.solstats = sol_stat;
            filename = strcat("./results/_ablation_8_", char(params_num.sim.aero_model), "_V_", num2str(baseWindSpeed), "_Nw_", num2str(params_num.sim.winding_number), ".mat");
            save(filename, "results");
        else
            results = struct();
            results.optimal_sol = optimal_resutls;
            results.time = t;
            results.solstats = sol_stat;
            filename = strcat("./results/_normal_8_", char(params_num.sim.aero_model), "_V_", num2str(baseWindSpeed), "_Nw_", num2str(params_num.sim.winding_number) , ".mat");
            save(filename, "results");
        end
    end
end