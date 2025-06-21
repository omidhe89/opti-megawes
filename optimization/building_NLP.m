% BUILDING_NLP.m
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

% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 14, 2025
% Last Modified: June 14, 2025
%
function [prob, prob_jacobians, lbw, ubw, lbg, ubg, w0, equality] = building_NLP(symbolic_funcs, x0, ocp_bounds, params_sym, params_num)
% Description:
%   This function constructs the Nonlinear Programming (NLP) problem for
%   the optimal control of an Airborne Wind Energy System (AWES). It
%   defines the decision variables, bounds, constraints, and the objective
%   function using CasADi for symbolic manipulation. The NLP is set up for
%   a multiple shooting formulation, incorporating state dynamics,
%   algebraic equations, and various path and control mode-specific constraints.
%
% Inputs:
%   - symbolic_funcs: A structure containing CasADi symbolic functions for
%     the system dynamics, algebraic constraints, and cost calculation.
%   - x0: Initial guess for the state, control, and algebraic variables,
%     and homotopy parameters.
%   - ocp_bounds: A structure containing lower and upper bounds for all
%     decision variables and system parameters in the Optimal Control Problem.
%   - params_sym: A structure containing symbolic representations of system parameters.
%   - params_num: A structure containing numerical values for system parameters.
%
% Outputs:
%   - prob: A CasADi structure defining the NLP problem, including the
%     objective function 'f', decision variables 'x', and constraints 'g'.
%   - prob_jacobians: A CasADi structure containing the symbolic Jacobians
%     of the objective function and constraints with respect to the
%     decision variables.
%   - lbw: Lower bounds for the decision variables.
%   - ubw: Upper bounds for the decision variables.
%   - lbg: Lower bounds for the constraints.
%   - ubg: Upper bounds for the constraints.
%   - w0: Initial guess for all decision variables.
%   - equality: A boolean array indicating which constraints are equality constraints.
%
% Dependencies:
%   - CasADi
%   - struct_to_key_value_array: A helper function to convert parameter
%     structures into a key-value array.
%   - scaled_to_si: A helper function to convert scaled variables to SI units.
%   - si_to_scaled: A helper function to convert SI units to scaled variables.
%   - model_for_optimization: A function that defines the system dynamics,
%     algebraic equations, and specific constraints for each time step.

import casadi.*


keys = fieldnames(params_sym);
params_num_array = [];
for i = 1:length(keys)
    params_num_field = struct_to_key_value_array(params_num.(keys{i}));
    params_num_array = [params_num_array; params_num_field];
end
% Convert each element of the second column to a column vector
params_column_vector = cellfun(@(x) x(:), params_num_array(:, 2), 'UniformOutput', false);

% Now concatenate them vertically (using vertcat)
paramsColumnVec = vertcat(params_column_vector{:});


Nx = params_num.Nx;    % x(1:2) -> winch , x(3:14) -> ap4
Nu = params_num.Nu;    % according to the control mode
Nz = params_num.Nz;    % [Magnitude, el, az];
N = params_num.sim.Nt;



xi_scaled = x0.xi_scaled;
ui_scaled = x0.ui_scaled;
zi_scaled = x0.zi_scaled;
phi_i = x0.phi_i;
T_i = x0.T_i;

phi_L = ocp_bounds.phi_L;
phi_U = ocp_bounds.phi_U;

x_L = ocp_bounds.x_L;
x_U = ocp_bounds.x_U;
z_L = ocp_bounds.z_L;
z_U = ocp_bounds.z_U;
u_L = ocp_bounds.u_L;
u_U = ocp_bounds.u_U;

% Start with an empty NLP

w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];
equality = [];


% homotopy parameters
if params_num.ctrl_mode == 3
    Phi = MX.sym('Phi', params_num.homotopy.phi_numbers);
    w = {w{:}, Phi};
    lbw = [lbw; phi_L];
    ubw = [ubw; phi_U];
    w0 = [w0;  phi_i];
else
    Phi = phi_i;
end



% "Lift" initial conditions
Xk = MX.sym('X_0', Nx);
w = {w{:}, Xk};
if params_num.sim.single_reel_out_activated == 1
    lbw = [lbw; [0; 0.5; zeros(Nx-2,1)]];
    ubw = [ubw; ones(Nx,1)];
else
    lbw = [lbw; zeros(Nx,1)];
    ubw = [ubw; ones(Nx,1)];
end
w0 = [w0;  xi_scaled(:,1)];


Zk = MX.sym('Z_0',Nz);
w = {w{:}, Zk};
lbw = [lbw;  zeros(Nz,1)];
ubw = [ubw;  ones(Nz,1)];
w0 = [w0;  zi_scaled(:,1)];

% cold start z0
Xk_si = scaled_to_si(Xk, x_L, x_U);
Zk_si = scaled_to_si(Zk, z_L, z_U);

alg_eq = symbolic_funcs.tether_alg_equation(Xk_si, Zk_si, paramsColumnVec);

g = {g{:}, vertcat(alg_eq)./ sqrt((z_U - z_L)/2)}; %
lbg = [lbg; zeros(Nz,1)];
ubg = [ubg; zeros(Nz,1)];
equality = [equality; true(Nz, 1)];
%
% constrain the init position
g = {g{:}, vertcat(Xk(12) - xi_scaled(12,1))}; 
lbg = [lbg; -0.0000 * ones(1,1)]; 
ubg = [ubg; 0.0000 * ones(1,1)]; 
equality = [equality; [true(1,1)]];


% Define time as a decision variable (to be optimized)
T = MX.sym('T');
% Create multiple shooting intervals
T_sim = (1 - Phi(3)) * scaled_to_si(T, ocp_bounds.time_L, ocp_bounds.time_U) + Phi(3) * params_num.sim.Tf;
h =  T_sim / N ; % Step size for each interval

% Formulate the NLP
X0 = Xk;
Xi = X0;
Zi = Zk;
J = 0;   % running cost
E = 0;   % mechanical energy E(T)/T -> terminal cost

for k=0:N-1

    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],Nu);
    w = {w{:}, Uk};
    lbw = [lbw; zeros(Nu,1)];
    ubw = [ubw; ones(Nu,1)];
    w0 = [w0;  ui_scaled(:,k+1)];


    Xk = MX.sym(['X_' num2str(k+1)],Nx);
    w = {w{:}, Xk};
    if params_num.sim.single_reel_out_activated == 1
        if  k <= floor(params_num.sim.reel_out_time*params_num.sim.Nt) % note >= Just here !
            lbw = [lbw; [0; 0.50; zeros(Nx-2,1)]];
            ubw = [ubw; ones(Nx,1)];
        else
            lbw = [lbw; zeros(Nx,1)];
            ubw = [ubw; ones(Nx,1)];
        end
    else
        lbw = [lbw; zeros(Nx,1)];
        ubw = [ubw; ones(Nx,1)];
    end
    w0 = [w0;  xi_scaled(:,k+2)];

    Zk = MX.sym(['Z_' num2str(k+1)],Nz);
    w = {w{:}, Zk};
    lbw = [lbw;  zeros(Nz,1)];
    ubw = [ubw;  ones(Nz,1)];
    w0 = [w0;  zi_scaled(:,k+2)];

    [Pk, omega_dot, lift_force, winch_acc, alpha, beta, g, lbg, ubg, equality] = model_for_optimization(symbolic_funcs, Xk, X0, Uk, Zk, Phi, g, lbg, ubg, params_num, paramsColumnVec, ocp_bounds, equality, h);

    % cost function
    beta_scaled = si_to_scaled(beta, ocp_bounds.beta_min, ocp_bounds.beta_max);
    alpha_scaled = si_to_scaled(alpha, ocp_bounds.alpha_min, ocp_bounds.alpha_max);
    J =  J +  h .* symbolic_funcs.cost(Xk, Zk, Uk, alpha_scaled, beta_scaled, omega_dot, lift_force, winch_acc, Phi, xi_scaled(:,k+2));
    E = E + h * Pk;

    % Ready for the next time step
    X0 = Xk;
end

% add mayer cost function (energy and time)
J = (J - 15 * (1 - Phi(3)) * E + 0.25 * (1 - Phi(3)) * (T_sim - params_num.sim.Tf)^2 )/ (N); % add terminal cost

% periodicty constraint
if params_num.sim.path_type == 8
    % To fly lemniscate
    g = {g{:}, vertcat(Xk - Xi)};
    lbg = [lbg; [zeros(14,1)]];
    ubg = [ubg; [zeros(14,1)]];
    equality = [equality; true(14, 1)];
else
    % To fly circule
    g = {g{:}, (vertcat(Xk(1:10) - Xi(1:10), Xk(11) - (si_to_scaled(scaled_to_si(Xi(11),x_L(11),x_U(11)) - 2*pi*params_num.sim.winding_number,x_L(11),x_U(11))) , Xk(12:end) - Xi(12:end)))}; % CCW
%      g = {g{:}, (vertcat(Xk(1:10) - Xi(1:10), Xk(11) - (si_to_scaled(scaled_to_si(Xi(11),x_L(11),x_U(11)) + 2*pi*params_num.sim.winding_number,x_L(11),x_U(11))) , Xk(12:end) - Xi(12:end)))}; %CW
    lbg = [lbg; [zeros(14,1)]];
    ubg = [ubg; [zeros(14,1)]];
    equality = [equality; [true(14, 1)]];
end

%
w = {w{:},  T};
lbw = [lbw; 0];
ubw = [ubw; 1];
w0 = [w0; si_to_scaled(T_i, ocp_bounds.time_L, ocp_bounds.time_U)];


% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
prob_jacobians = struct('f', jacobian(J,vertcat(w{:})), 'g', jacobian(vertcat(g{:}),vertcat(w{:})));
end