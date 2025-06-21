function T_final = guess_final_time(params_num, init_path_radius)
% GUESS_FINAL_TIME Estimates the final simulation time for MegAWES OCP.
%
% Copyright (c) 2025 [Omid Heydarnia]
% This function is part of the 'Opti-MegAWES' tool and is licensed under the
% Apache License, Version 2.0. See the accompanying LICENSE file for the
% full text of the license.
%
% Author: [Omid Heydarnia]
% Date: June 12, 2025
% Last Modified: June 12, 2025
%
% Description:
%   This function provides an educated initial guess for the total simulation
%   time (T_final) required for the MegAWES (Mega Airborne Wind Energy System)
%   Optimal Control Problem (OCP). The estimation is based on the selected
%   flight path type (circle or lemniscate), the desired number of winding
%   cycles, and an initial estimate of the flight speed and path radius.
%   It calculates the approximate total path length (circumference) and
%   divides it by the initial flight speed.
%
% Inputs:
%   params_num       - Structure containing numerical parameters for the
%                      simulation, specifically:
%                      .sim.path_type      (0 for circle, 8 for lemniscate)
%                      .sim.init_flight_speed (m/s, initial guess for kite speed)
%                      .sim.winding_number (Number of flight cycles)
%   init_path_radius - double, Initial radius (m) used for the flight path
%                      (e.g., radius of circle, characteristic radius of lemniscate).
%
% Outputs:
%   T_final          - double, Estimated final time (s) for the simulation,
%                      rounded down to the nearest integer.
%
% Example:
%   params.sim.path_type = 0;
%   params.sim.init_flight_speed = 90;
%   params.sim.winding_number = 3;
%   radius = 240;
%   T_sim = guess_final_time(params, radius);
%
% See also: SOLVE_MEGAWES_OCP

% Extract relevant numerical parameters from the input structure for clarity.
path_type = params_num.sim.path_type;
init_flight_speed = params_num.sim.init_flight_speed; 
winding_number = params_num.sim.winding_number; 

if path_type == 0 
    circumference = 2 * pi * winding_number * init_path_radius;
else
    circumference = 4 * 0.9 * pi * winding_number * init_path_radius; % roughly estimates the lemniscate with two circle
end
T_final = floor(circumference / init_flight_speed);
T_final = floor(T_final); %0.05 
end