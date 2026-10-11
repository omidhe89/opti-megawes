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

% This code is an interpretation and adaptation of concepts from objective.py
% from awebox toolbox developed in freiburg university. This file is part of awebox.
% awebox -- A modeling and optimization framework for multi-kite AWE systems.

function cost = cost_function(vars, params_num)
% COST_FUNCTION Defines the objective function for the optimal control problem.
%
% This function constructs the total cost (objective) to be minimized in the
% optimal control problem for a kite power system. The cost function is a
% weighted sum of several terms, including:
%   - **Homotopy Cost**: Used in homotopy methods to guide the solver.
%   - **Initialization Cost**: Penalizes winch and angular acceleration to
%     avoid aggressive manoeuvre.
%     initial accelerations, typically active at the first collocation point.
%   - **Tracking Cost**: Penalizes deviations from a initial reference trajectory
%     and potentially encourages specific aerodynamic conditions (e.g., lift).
%   - **Power Cost**: Maximizes generated power.
%
% The terms are weighted by homotopy parameters (phi) to allow for a gradual
% introduction of different objectives during the optimization process.
%
% Inputs:
%   vars            - Struct, containing symbolic CasADi variables:
%                     .x              - State vector (scaled).
%                     .u              - Control input vector (scaled).
%                     .z              - Algebraic variable vector (scaled).
%                     .phi            - Homotopy parameter vector.
%   params_num      - Struct, numerical parameters.
%
% Outputs:
%   cost            - CasADi Function object, representing the total cost.
%                     It takes unscaled state (x), algebraic (z), control (u),
%                     scaled alpha/beta, scaled omega_dot, unscaled lift_force,
%                     scaled winch_acc, homotopy parameters (phi), and scaled
%                     reference state (xi_scaled) as inputs.
%


% Import CasADi library for symbolic computations.
import casadi.*

% Declare symbolic variables that are not part of the 'vars' struct but are
% inputs to the cost function for various penalty terms.
% These are typically scaled versions or derived quantities.
beta_scaled = casadi.SX.sym('beta_scaled');     % Scaled sideslip angle.
alpha_scaled = casadi.SX.sym('alpha_scaled');   % Scaled angle of attack.
xi_scaled = casadi.SX.sym('xi_scaled', params_num.Nx); % Scaled reference state vector for tracking.
omega_dot = casadi.SX.sym('omega_dot', 3);     % Angular acceleration vector (scaled) for initial cost.
v_dot = casadi.SX.sym('v_dot', 3);     % linear acceleration vector (scaled) for initial cost.
lift_force = casadi.SX.sym('lift_force', 1);   % Lift force vector (unscaled) for tracking cost.
winch_acc = casadi.SX.sym('winch_acc', 1);     % Winch acceleration (scaled) for initial cost.


% Extract actual (unscaled) variables from the input 'vars' struct.
x = vars.x; % State vector
u = vars.u; % Control input vector
z = vars.z; % Algebraic variable vector
phi = vars.phi; % Homotopy parameters

% Scale the lift force for inclusion in the cost function.
% The scaling factor (1e5) normalizes the lift force to a comparable magnitude
% with other terms, preventing it from dominating or being negligible.
lift_force_scaled = lift_force ./ 1e5;

% Calculate the error between the current state and the reference state.
err = x - xi_scaled;

% --- Homotopy Cost Term ---

if params_num.ctrl_mode == 3
    % In control mode 3, a linear penalty on the sum of homotopy parameters.
    % The coefficient 0.08 'turns on' this cost as phi increases.
    cost_homotopy = 0.08 * ones(1, params_num.homotopy.phi_numbers) * phi;
    % --- Tracking Cost Term ---
    % weights x1..x15; x12:14 = W-frame position (x,y,z) and x15 = q0
    cost_tracking = err' * 1.05 * diag([0, 0.000, 1.5, 0.000 * ones(1, 2), 0.000 * ones(1, 3), 0.000 * ones(1, 3), 0.61551, 0.010551, 0.01505, 0.010]) * err; %- ...
                %(1 - phi(1)) * 3.5 * (lift_force_scaled' * lift_force_scaled);
else
    % For other control modes, the homotopy and tracking costs are zero, implying
    % no homotopy continuation is being applied to guide the optimization.
    cost_homotopy = 0;
    cost_tracking = 0;
end

% --- Initialization Cost Term ---
cost_init = omega_dot' * 25.5 * diag([1.3525, 1.6525, 1.5525]) * omega_dot + ...
            v_dot' * 7.5 * diag([1.9525, 0.6525, 0.6525]) * v_dot + ...
            5.500 * beta_scaled^2 ;% + 0.5000 * winch_acc^2 - 0.8* x(1)^2; + ...

% --- Power Cost Term ---
% This term aims to maximize the generated power.
P = params_num.winch.radius * z(1) * x(2);
cost_power = -40 * P;

% --- Total Objective Function (J) ---
J = cost_init + phi(2) * cost_tracking + phi(3) * cost_power + cost_homotopy;

% Create a CasADi Function object for the cost.
% This allows the cost function to be efficiently evaluated within the
% optimization solver. The inputs and outputs are explicitly defined.
cost = casadi.Function('cost', ...
                       {x, z, u, alpha_scaled, beta_scaled, omega_dot, v_dot, lift_force, winch_acc, phi, xi_scaled}, ...
                       {J}, ...
                       {'x_scaled', 'z_scaled', 'u_scaled', 'alpha_scaled', 'beta_scaled', 'omega_dot_scaled', 'v_dot_scaled', 'lift_force_unscaled', 'winch_acc_scaled', 'phi_parameters', 'xi_scaled_reference'}, ...
                       {'J'}); % Output name for the objective function.

end