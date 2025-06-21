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
% Email: omid.heydarnia@hotmail.com
% Date: June 12, 2025
% Last Modified: June 12, 2025

function w0 = NLP_initialization(init_mode, symbolic_funcs, params_sym, params_num, ocp_bounds, init_path_radius)
    % NLP_INITIALIZATION Initializes the NLP variables for optimal control problem.
    %
    % This function sets up the initial guess for the state (xi), algebraic
    % variables (zi), and control inputs (ui) for the Nonlinear Programming (NLP)
    % problem. It supports different initialization modes (single-loop or multi-loop)
    % and control modes, and also scales the initial guess according to bounds.
    %
    % Inputs:
    %   init_mode         - String, Specifies the initialization strategy:
    %                       'single' for single-loop path initialization.
    %                       'multi' for multi-loop path initialization (via path_init_mlt_loop, still in progress).
    %   symbolic_funcs    - Struct, Containing symbolic functions required for path initialization.
    %   params_sym        - Struct, Symbolic parameters of the system.
    %   params_num        - Struct, Numerical parameters of the system, 
    %   ocp_bounds        - Struct, Contains lower (x_L, z_L, u_L) and upper (x_U, z_U, u_U)
    %                       bounds for states, algebraic variables, and controls.
    %   init_path_radius  - Scalar, Initial path radius for circular trajectory initialization.
    %
    % Outputs:
    %   w0                - Struct, Containing the initial guess for NLP variables:
    %                       w0.xi_scaled    - Initial guess for states, scaled between 0 and 1.
    %                       w0.zi_scaled    - Initial guess for algebraic variables, scaled.
    %                       w0.ui_scaled    - Initial guess for control inputs, scaled.
    %                       w0.phi_i        - Initial guess for homotopy parameters.
    %                       w0.T_i          - Initial guess for final time.

    w0 = struct();

    % Determine the path initialization strategy based on 'init_mode'.
    if strcmp(init_mode, 'single')
        % Initialize path for a single loop.
        [xi, zi, ui] = path_init(symbolic_funcs, params_sym, params_num, init_path_radius);
    else % Assumed 'multi' or other mode for multi-loop initialization
        % Initialize path for multiple loops.
        [xi, zi, ui] = path_init_mlt_loops(symbolic_funcs, params_sym, params_num, init_path_radius); % ( still in progress!)
    end

    % Adjust the initial control input (ui) based on the control mode.
    if params_num.ctrl_mode == 1
        ui = [ui(1,:); ui(2:7,:); zeros(3, params_num.sim.Nt); zeros(1, params_num.sim.Nt)];
    elseif params_num.ctrl_mode == 2
        ui = [ui(1,:); ui(2:end,:)];
    else
        ui = [ui(1,:); ui(2:7,:); zeros(3, params_num.sim.Nt); ocp_bounds.u_U(11) .* ones(1, params_num.sim.Nt); ui(2:7,:)];
    end

    % Initialize homotopy parameters (phi_i).
    phi_i = ones(params_num.homotopy.phi_numbers, 1);

    % Scale the initial conditions for optimization.
    N = params_num.sim.Nt; % Number of time intervals
    for i = 0:N
        % Scale states (xi)
        w0.xi_scaled(:, i + 1) = si_to_scaled(xi(:, i + 1), ocp_bounds.x_L, ocp_bounds.x_U);
        % Scale algebraic variables (zi)
        w0.zi_scaled(:, i + 1) = si_to_scaled(zi(:, i + 1), ocp_bounds.z_L, ocp_bounds.z_U);
        if i > 0
            % Scale control inputs (ui) (controls are defined at interval starts, hence i > 0)
            w0.ui_scaled(:, i) = si_to_scaled(ui(:, i), ocp_bounds.u_L, ocp_bounds.u_U);
        end
    end

    % Set initial values for homotopy parameters based on control mode.
    if params_num.ctrl_mode == 1
        % If direct control, homotopy parameters are often initialized to zero,
        % implying a direct solution without a homotopy path.
        phi_i = zeros(params_num.homotopy.phi_numbers, 1);
    else
        phi_i = ones(params_num.homotopy.phi_numbers, 1);
    end
    w0.phi_i = phi_i;

    % Set the initial guess for the final time (T_i).
    w0.T_i = params_num.sim.Tf;
end