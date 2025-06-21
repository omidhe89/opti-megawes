% Copyright (c) 2025 Omid Heydarnia
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 14, 2025
% Last Modified: June 14, 2025

function latex_description(latexFile, params_num, ocp_bounds)
% to write problem formulation, variable bounds, and OCP constraints.
%   latexFile: the latex file we write in
%   params_num: system and flight parameters
%   ocp_bounds: the variable bounds and constraints limits



fprintf(latexFile, 'We aim to generate periodic maximum power flight trajectories for an Airborne Wind Energy System (AWES). Similar to previous work, we propose to solve the following periodic continuous-time Optimal Control Problem (OCP): \n');

fprintf(latexFile, '\\begin{equation}\n');
fprintf(latexFile, '\\label{eq:ocp}\n');
fprintf(latexFile, '\\begin{aligned}\n');
fprintf(latexFile, '&\\min_{\\vec{x}, \\vec{z}, \\vec{u}, T} \\ \\frac{1}{T} \\int_{0}^{T} c(\\vec{x}(t), \\vec{u}(t), \\vec{z}(t)) \\, \\mathrm{d}t \\\\\n');
fprintf(latexFile, '\\text{s.t.} \\quad \n');
fprintf(latexFile, '&0 = \\vec{F}(\\dot{\\vec{x}}(t), \\vec{x}(t), \\vec{u}(t), \\vec{z}(t)), \\quad \\forall t \\in [0, T) \\\\\n');
fprintf(latexFile, '&0 \\geq \\vec{H}(\\vec{x}(t), \\vec{u}(t), \\vec{z}(t)), \\quad \\forall t \\in [0, T) \\\\\n');
fprintf(latexFile, '&0 = \\vec{x}(0) - \\vec{x}(T) \\\\\n');
fprintf(latexFile, '&0 = \\Psi(\\vec{x}(0))\n');
fprintf(latexFile, '\\end{aligned}\n');
fprintf(latexFile, '\\end{equation}\n');

fprintf(latexFile, 'The OCP is fully specified by the definition of the cost rate function $c$, the differential and other equality constraints $\\vec{F}$, and the inequality path constraints $\\vec{H}$.\n');
fprintf(latexFile, '\n Some of important variable bounds and path constraints are introduced in the followig table \n');

table_format_start = [ ...
'\\begin{table}[ht!]\n' ...
'\\centering\n' ...
'\\caption{variable bounds and path constraints}\n' ...
'\\begin{tabular}{ c c c c c }\n' ...
'\\hline\n' ...
'Description & Symbol & Min & Max & Unit \\\\ \n' ...
'\\hline\n'];
fprintf(latexFile, table_format_start);
table_raw = sprintf('Apparent velocity & $\\|v_\\text{a}\\|$ & $%.1f$  & $%.1f$  & $m/s$ \\\\ \n', ocp_bounds.v_app_min, ocp_bounds.v_app_max);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('Angle of attack & $\\alpha$ & $%.1f$  & $%.1f$  & $deg$ \\\\ \n', ocp_bounds.alpha_min *180/pi, ocp_bounds.alpha_max * 180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('Side slip angle & $\\beta$ & $%.1f$  & $%.1f$  & $deg$ \\\\ \n', ocp_bounds.beta_min *180/pi, ocp_bounds.beta_max * 180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('Tether force (aircraft) & $\\|F_{\\text{tether}}\\|$ & $%.2f\\times10^3$  & $%.3f\\times10^6$  & $N$ \\\\ \n', ocp_bounds.z_L(1)/1e3, ocp_bounds.z_U(1)/1e6);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('Winch acceleration & $\\dot{\\omega_w}$ & $%.1f$  & $%.1f$  & $rad/s^2$ \\\\ \n', ocp_bounds.winch_acc_min, ocp_bounds.winch_acc_max);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('aileron deflection & $\\delta_{ai}$ & $%.1f$  & $%.1f$  & $deg$ \\\\ \n', ocp_bounds.u_L(8)*180/pi, ocp_bounds.u_U(8)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('elevator deflection & $\\delta_{el}$ & $%.1f$  & $%.1f$  & $deg$ \\\\ \n', ocp_bounds.u_L(9)*180/pi, ocp_bounds.u_U(9)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('rudder deflection & $\\delta_{ru}$ & $%.1f$  & $%.1f$  & $deg$ \\\\ \n', ocp_bounds.u_L(10)*180/pi, ocp_bounds.u_U(10)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('angular velocity & $\\omega_{x}$ & $%.1f$  & $%.1f$  & $deg/s$ \\\\ \n', ocp_bounds.x_L(6)*180/pi, ocp_bounds.x_U(6)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('angular velocity & $\\omega_{y}$ & $%.1f$  & $%.1f$  & $deg/s$ \\\\ \n', ocp_bounds.x_L(7)*180/pi, ocp_bounds.x_U(7)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_raw = sprintf('angular velocity & $\\omega_{z}$ & $%.1f$  & $%.1f$  & $deg/s$ \\\\ \n', ocp_bounds.x_L(8)*180/pi, ocp_bounds.x_U(8)*180/pi);
fprintf(latexFile, '%s', table_raw);
table_format_end = [ ...
'\\hline\n' ...
'\\end{tabular}\n' ...
'\\label{tab:vars_constraints}\n' ...
'\\end{table}\n' ...
];
fprintf(latexFile, table_format_end);

fprintf(latexFile, '\n An optimal trajectory with \\textbf{%d} pumping cycle is generated for MegAWES aircraft where the \\textbf{%s}  aerodynamic model is utilized. The wind speed is \\textbf{%d} $m/s$.', params_num.sim.winding_number , params_num.sim.aero_model, params_num.env.wind.base_velocity);
end