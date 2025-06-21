% SOLVE_MEGAWES_OCP.m
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% See the accompanying LICENSE file for the full text of the license.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 10, 2025
% Last Modified: June 10, 2025
%
% 
% The implementation of the penalty-based homotopy approach in this script is
% interpreted and adapted from the methodology described in "De Schutter, J.;
% Leuthold, R.; Bronnenmeyer, T.; Malz, E.; Gros, S.; Diehl, M. AWEbox: An Optimal
% Control Framework for Single- and Multi-Aircraft Airborne Wind Energy Systems.
% Energies 2023, 16, 1900. https://doi.org/10.3390/en16041900". This approach
% facilitates the solution of complex optimal control problems by gradually
% transforming an easier-to-solve problem into the desired, more challenging one
% through a series of intermediate optimizations, often improving convergence.
%
% Description:
%   This program is designed to solve an Optimal Control Problem (OCP) for a
%   MegAWES (Mega Airborne Wind Energy System) model. It uses CasADi for
%   symbolic differentiation and IPOPT for numerical optimization to find
%   optimal flight trajectories and control inputs for power generation.
%   The script supports different tether models (flexible/rigid), path types
%   (circle/lemniscate), and control modes (aerodynamic, fictitious, homotopy).
%
% Inputs:
%   None directly; parameters are set within the script.
%
% Outputs:
%   - Optimal states (position, velocity, etc.) and control inputs for the MegAWES.
%   - Plots and animations of the optimal flight trajectory.
%   - Saved flight data for further analysis.
%
% Dependencies:
%   - CasADi (must be installed and on MATLAB path)
%   - Custom libraries: 'lib', 'MegAWES', 'others', 'transformations',
%     'environment', 'optimization', 'results' (these directories must
%     contain necessary functions for model definition, OCP formulation,
%     initialization, and visualization).

%% Initialization and Path Setup
% Clear workspace, command window, and add necessary directories to MATLAB path.
close all
clear
clc
addpath('lib','MegAWES','others','transformations','environment', 'optimization','results')

% Import CasADi library for symbolic computation.
import casadi.*

%% Set Simulation and Model Parameters
% Define key numerical parameters for the MegAWES model and the OCP simulation.
baseWindSpeed = 15;     % m/s, reference wind speed at a specific altitude. 
baseWindDir = pi;       % radians, wind direction. Pi (180 degrees) means
                        % wind flows along the positive x-axis in a standard Cartesian coordinate system.

% Simulation and Path Parameters
% Note: The circle flies counter-clockwise
params_num.sim.ablation_study = 0;      % 0: flexible tether, 1: rigid tether
path_type = 0;          % Defines the desired flight path geometry.
                        % 0: circular path, 8: lemniscate (figure-eight) path.
init_path_radius = 220; % m, initial radius for path generation (for both circle and lemniscate).
params_num.sim.path_type = path_type;
params_num.sim.init_flight_speed = 85; % m/s, initial guess for the aircraft's flight speed.
params_num.sim.winding_number = 3; % cylces numbers ( Note: winding number for circles should be twice the size of lemniscates for comparison.) 

% Calculate Estimated Final Time (Tf) for the flight trajectory.
% Tf estimates an appropriate simulation duration based on path type and winding number.
params_num.sim.Tf = guess_final_time(params_num, init_path_radius);
msg = sprintf("the proposed time for flight is: %0.5g", params_num.sim.Tf);
disp(msg);

% Determine the number of shooting nodes (discretization points) for the OCP.
% This impacts the resolution and computational cost of the optimization.
estimated_time_step = 0.5;
params_num.sim.Nt = floor(params_num.sim.Tf / estimated_time_step) + floor(0.5 * params_num.sim.Tf / estimated_time_step) + 1;
msg = sprintf("the proposed shooting nodes for flight is: %0.5g", params_num.sim.Nt);
disp(msg);

params_num.sim.single_reel_out_activated = 1; % Reel-out strategy:
                                              % 0: multiple reel-out phases, 1: single reel-out phase.
params_num.sim.reel_out_time = 0.7; % Normalized time (0 to 1) at which reel-out occurs if single_reel_out_activated is 1.

params_num.sim.aero_model = 'ALM';  % Aerodynamic model choice for the kite:
                                    % 'ALM' (Actuator Line Method),
                                    % 'VLM' (Vortex Lattice Method),
                                    % 'CFD' (Computational Fluid Dynamics),
                                    % 'MegAWES' (specific MegAWES model).

% Control Mode:
% 1: Aerodynamic control (direct manipulation of aerodynamic forces/moments).
% 2: Fictitious control (abstract control inputs without direct physical mapping).
% 3: Homotopy control (uses a continuation method to gradually solve the OCP).                                    
params_num.ctrl_mode = 3;
% Homotopy parameters (if 'ctrl_mode' is 3).
params_num.homotopy.phi_numbers = 3; % Number of homotopy steps for progressive optimization. (always 3)

%% Symbolic Model Generation
% Begin the process of building the symbolic representation of the MegAWES
% system dynamics and the optimal control problem. This step uses CasADi
% to define variables, parameters, and system equations.
msg = newline + "start building symbolic model ...";
disp(msg);
tic
% define symbolic variables (states, controls, algebraic variables) and symbolic/numeric parameters of the system.
[vars, params_sym, params_num] = system_vars_params(baseWindSpeed, baseWindDir, params_num);

%% Set OCP Variables Bounds
% Define the lower and upper bounds for the decision variables (states, controls, and algebraic variables) in the Optimal Control Problem.
ocp_bounds = system_OCP_bounds(params_num);

%% Build Symbolic System Dynamics and Cost Function
% Define the symbolic equations governing the MegAWES system and the objective (cost) function to be minimized/maximized.
symbolic_funcs = struct();
symbolic_funcs.data.params_sym = params_sym; 
% components of the MegAWES model:
% - kite_data: functions related to kite position and orientation.
% - tether_dynamics: tether quantities (e.g. forces at winch and kite, tether lumped mass)
% - tether_alg_equation: algebraic equations related to tether constraints.
% - kite_dynamics: differential equations for kite motion.
% - winch_dynamics: differential equations for winch system.
% - aerodynamics_vars: aerodynamic quantities (e.g., airspeed, angles).
% - F_M_aero: aerodynamic forces and moments.
% - data.array_params_sym: symbolic parameters arranged in an array format.
[symbolic_funcs.kite_data, symbolic_funcs.tether_dynamics, symbolic_funcs.tether_alg_equation, ...
    symbolic_funcs.kite_dynamics, symbolic_funcs.winch_dynamics, symbolic_funcs.aerodynamics_vars,...
    symbolic_funcs.F_M_aero, symbolic_funcs.data.array_params_sym] = build_dynamics(vars, params_sym, params_num);
toc

%% OCP Initialization
% Create an initial guess for circular and lemniscate trajectories.
% Finding the initial states (differential and algebraic) and control inputs for optimization.
msg = newline + "start initialization ....";
disp(msg);
tic
w_init = NLP_initialization('single', symbolic_funcs, params_sym, params_num, ocp_bounds, init_path_radius); % choose between 'single' and 'multi' ( multi is in progress currently)
toc
%% Building the Nonlinear Program (NLP)
% Formulate the Optimal Control Problem as a Nonlinear Program. 
msg = newline + "start building NLP ....";
disp(msg)
tic

% Define the lagrange cost function to be evaluated through the entire time horizon.
symbolic_funcs.cost = cost_function(vars, params_num);
% Construct the CasADi NLP problem including:
% prob: CasADi NLP, containing the objective and constraints.
% prob_jacobians: Structure containing symbolic Jacobians (if needed).
% lbw, ubw: Lower and upper bounds for the decision variables (w).
% lbg, ubg: Lower and upper bounds for the constraints (g).
% w0: Initial guess for the decision variables.
% equality: Structure to check equality and in-equality constraints.
[prob, prob_jacobians, lbw, ubw, lbg, ubg, w0, equality] = building_NLP(symbolic_funcs, w_init, ocp_bounds, params_sym, params_num);
 
% Configure IPOPT solver options based on the chosen control mode.
if params_num.ctrl_mode == 3
    % Configure optimization parameters for the initial (most trivial) OCP.
    % This solver aims for quick convergence to a feasible, even if not optimal, solution.
    opts_init = struct();
    opts_init.ipopt.max_iter = 3000;
    opts_init.ipopt.mu_init = 1;   % tau_0
    opts_init.ipopt.mu_target = 0.01; % tau_i
    opts_init.ipopt.tol = 1e-8;
    opts_init.ipopt.acceptable_iter = 5;
    opts_init.ipopt.sb = char('yes'); % suppress IPOPT banner
    opts_init.ipopt.print_level = 5;
    opts_init.ipopt.nlp_scaling_method = char('none');
    opts_init.expand = true;
    solver_init = nlpsol('solver', 'ipopt', prob, opts_init);
    % inermediate homotopy solvers
    % slightly looser tolerance for faster intermediate convergence.
    opts_warmstart = struct();
    opts_warmstart.ipopt.max_iter = 500;
    opts_warmstart.ipopt.mu_init = 0.01;   % tau_0
    opts_warmstart.ipopt.mu_target = 0.01; % tau_i
    opts_warmstart.ipopt.tol = 1e-5;
    opts_warmstart.ipopt.acceptable_iter = 15;
    opts_warmstart.ipopt.sb = char('yes'); % suppress IPOPT banner
    opts_warmstart.ipopt.print_level = 5;
    opts_warmstart.ipopt.warm_start_init_point = char('yes');
    opts_warmstart.ipopt.nlp_scaling_method = char('none');
    opts_warmstart.expand = true;
    solver_warmstart = nlpsol('solver', 'ipopt', prob, opts_warmstart);
    % final solver
    % Options for the final, most refined homotopy step.
    % Aims for high accuracy and tight convergence.
    opts_final = struct();
    opts_final.ipopt.max_iter = 3000;
    opts_final.ipopt.mu_init = 0.01;
    opts_final.ipopt.mu_target = 0;
    opts_final.ipopt.tol = 1e-8;
    opts_final.ipopt.constr_viol_tol = 1e-8;
    opts_final.ipopt.sb = char('yes'); % suppress IPOPT banner
    opts_final.ipopt.print_level = 5;
    opts_final.ipopt.warm_start_init_point = char('yes');
    opts_final.ipopt.nlp_scaling_method = char('none');
    opts_final.expand = true;
    solver_final = nlpsol('solver', 'ipopt', prob, opts_final);
    
    toc
    
    msg = newline + " start optimization ....";
    disp(msg)
    %% Solve the NLP using Homotopy Method
    % The homotopy method progressively introduces complexity into the OCP,
    % starting with a simpler problem and gradually transitioning to the full problem, helping the solver converge.
    tic
    disp('********  initial solver **********')
    sol = solver_init('x0', w0, 'lbx', lbw, 'ubx', ubw,...
                'lbg', lbg, 'ubg', ubg);
    T = toc;
    sol_stat = [];
    w_opt_0 = full(sol.x);
    stats = solver_init.stats();
    inf_pr_init = stats.iterations.inf_pr(end);
    inf_du_init = stats.iterations.inf_du(end);
    sol_stat = [sol_stat;T inf_pr_init inf_du_init];
    w_i = w_opt_0;
    for i=1:params_num.homotopy.phi_numbers
        lbw(i) = 0;
        tic
        formatSpec = "******** lower intermediate solver : %d, %s";
        txt = sprintf(formatSpec , i,  "********" )
        sol = solver_warmstart('x0', w_i, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg);
        T = toc;
        w_i = full(sol.x);
        stats = solver_warmstart.stats();
        inf_pr = stats.iterations.inf_pr(end);
        inf_du = stats.iterations.inf_du(end);
        sol_stat = [sol_stat;T inf_pr inf_du];
        ubw(i) = 0;
        
        tic
        formatSpec = "******** upper intermediate solver : %d, %s";
        txt = sprintf(formatSpec , i,  "********" )
        sol = solver_warmstart('x0', w_i, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg);
        T=toc;
        w_i = full(sol.x);
        stats = solver_warmstart.stats();
        inf_pr = stats.iterations.inf_pr(end);
        inf_du = stats.iterations.inf_du(end);
        sol_stat = [sol_stat;T inf_pr inf_du];
    end
    % Solve the final NLP with the 'solver_final'.
    tic
    disp('********  final solver **********')
    sol = solver_final('x0', w_i, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
    T = toc;
    w_opt = full(sol.x);
    stats = solver_final.stats();
    inf_pr_final = stats.iterations.inf_pr(end);
    inf_du_final = stats.iterations.inf_du(end);
    sol_stat = [sol_stat; T inf_pr_final inf_du_final];

else % Other Control Modes (e.g., Aerodynamic, Fictitious)
    %% Solve the NLP for non-homotopy modes
    sol_stat = [];
    
    opts_final = struct();
    opts_final.ipopt.max_iter = 5000;
    opts_final.ipopt.mu_init = 1;
    opts_final.ipopt.mu_target = 0;
    opts_final.ipopt.tol = 1e-9;
    opts_final.ipopt.constr_viol_tol = 1e-8;
    opts_final.ipopt.sb = char('yes'); % suppress IPOPT banner
    opts_final.ipopt.print_level = 5;
%     opts_final.ipopt.bound_relax_factor = 1e-8;
    opts_final.ipopt.warm_start_init_point = char('yes');
    opts_final.ipopt.nlp_scaling_method = char('none');
    opts_final.expand = true;
    
    solver_final = nlpsol('solver', 'ipopt', prob, opts_final);
    toc   
        
    msg = newline + " start optimization ....";
    disp(msg)
    tic
    sol = solver_final('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
    T = toc;
    w_opt = full(sol.x);
    stats = solver_final.stats();
    inf_pr_final = stats.iterations.inf_pr(end);
    inf_du_final = stats.iterations.inf_du(end);
    sol_stat = [sol_stat; T inf_pr_final inf_du_final];
    toc
end
%% Numerical Manipulation and Extracting Results
% Extract and scale to si the optimized states, algebraic variables, and control
% inputs from the raw 'w_opt' vector into meaningful arrays.
Nx = params_num.Nx;
Nu = params_num.Nu;
Nz = params_num.Nz;
N = params_num.sim.Nt; 
if (params_num.ctrl_mode == 1 || params_num.ctrl_mode == 2)
    un_scaled = zeros(Nu,N);
    zn_scaled = zeros(Nz,N+1);
    xn_scaled = zeros(Nx,N+1);
    total_len = Nx + Nu + Nz;
    for i=0:N   
        xn_scaled(:,i+1) = w_opt(total_len*i+1:total_len*i + Nx,1);
        zn_scaled(:,i+1) = w_opt(total_len*i+(Nx+1):total_len*i+(Nx+Nz),1);
        if i < N
            un_scaled(:, i+1) = w_opt(total_len*i+(Nx+Nz+1):total_len*i+(Nx+Nz+Nu),1);  
        end
        
    end
else
    un_scaled = zeros(Nu,N);
    zn_scaled = zeros(Nz,N+1);
    xn_scaled = zeros(Nx,N+1);
    phi_n = w_opt(1:params_num.homotopy.phi_numbers,1);
    total_len = Nx + Nu + Nz;
    for i=0:N   
        xn_scaled(:,i+1) = w_opt(total_len*i+(params_num.homotopy.phi_numbers+1):total_len*i+(Nx+params_num.homotopy.phi_numbers),1);
        zn_scaled(:,i+1) = w_opt(total_len*i+(Nx+params_num.homotopy.phi_numbers+1):total_len*i+(Nx+Nz+params_num.homotopy.phi_numbers),1);
        if i < N
            un_scaled(:, i+1) = w_opt(total_len*i+(Nx+Nz+params_num.homotopy.phi_numbers+1):total_len*i+(Nx+Nz+Nu+params_num.homotopy.phi_numbers),1);
        end
    end
end


t = (0:1/N:1).*  scaled_to_si(w_opt(end), ocp_bounds.time_L, ocp_bounds.time_U);
h = scaled_to_si(w_opt(end), ocp_bounds.time_L, ocp_bounds.time_U)/params_num.sim.Nt;
xn = scaled_to_si(xn_scaled, ocp_bounds.x_L, ocp_bounds.x_U);
un = scaled_to_si(un_scaled, ocp_bounds.u_L, ocp_bounds.u_U);
zn = scaled_to_si(zn_scaled, ocp_bounds.z_L, ocp_bounds.z_U);

pos_O = casadi.DM.zeros(3,N+1);
pos_W = casadi.DM.zeros(3,N+1);
vel_W = casadi.DM.zeros(3,N+1);
alpha = casadi.DM.zeros(1,N+1);
beta = casadi.DM.zeros(1,N+1);
M_OB = casadi.DM.eye(3);
% Prepare numerical parameters in an array format for symbolic functions.
% 'struct_to_key_value_array' converts nested structures to a flat array.
keys = fieldnames(params_sym);
params_num_array = [];
for i = 1:length(keys)
    params_num_field = struct_to_key_value_array(params_num.(keys{i}));
    params_num_array = [params_num_array; params_num_field];
end
% Convert each element of the second column to a column vector required for evaluating the CasADi symbolic functions.
params_column_vector = cellfun(@(x) x(:), params_num_array(:, 2), 'UniformOutput', false);
% Concatenate them vertically (using vertcat)
paramsColumnVec = vertcat(params_column_vector{:});
optimal_results = struct;
% Loop through each time step to compute additional derived quantities.
for k = 1:N+1
        optimal_results(k).xn = xn(:,k);
        optimal_results(k).zn = zn(:,k);
        [pos_O(:,k), pos_W(:,k), vel_W(:,k), M_OB, ~] = symbolic_funcs.kite_data(xn(:,k), paramsColumnVec);
        optimal_results(k).pos_O = full(pos_O(:,k));
        optimal_results(k).pos_W = full(pos_W(:,k));
        optimal_results(k).vel_W = full(vel_W(:,k));
        [Va, alpha(k), beta(k)] = symbolic_funcs.aerodynamics_vars(xn(:,k), paramsColumnVec);
        optimal_results(k).Va = Va;
        optimal_results(k).alpha = alpha(k);
        optimal_results(k).beta = beta(k);
        [err, T_kite, tether_pos] = symbolic_funcs.tether_dynamics(xn(:,k), zn(:,k), paramsColumnVec);
        optimal_results(k).tether_pos = full(tether_pos);
        p_opt = params_num.winch.radius*zn(1,k)*xn(2,k);
        optimal_results(k).p_opt = p_opt;
        T_kite_B = M_OB'*transformFromWtoO(params_num.env.wind.direction, full(T_kite));
        optimal_results(k).T_kite_B = T_kite_B;
        if k < N+1
            optimal_results(k).un = un(:,k);
        end
end


%% Visualization and Data Saving
% Generate plots and animations of the optimal flight trajectory and save
% relevant data for future analysis.
ui = scaled_to_si(w_init.ui_scaled, ocp_bounds.u_L, ocp_bounds.u_U);
xi = scaled_to_si(w_init.xi_scaled, ocp_bounds.x_L, ocp_bounds.x_U);
zi = scaled_to_si(w_init.zi_scaled, ocp_bounds.z_L, ocp_bounds.z_U);

% save optimal flight data in mat files
save_flight_data(t, optimal_results, sol_stat, params_num)
% create plots and save them in pdf format in the results directory
awe_plots(t, optimal_results, sol_stat, h, params_num, ocp_bounds)
% show animation of kite flying optimal trajectory
awe_animation(t, optimal_results, h, params_num);


