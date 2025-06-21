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
function [xi, zi, ui] = path_init(symbolic_funcs, params_sym, params_num, circRadius)
    % PATH_INIT Initializes a single-loop reference trajectory for the kite system.
    % For multiple loops trajectoris, kite flies a single loop multiple times. 
    % This function generates an initial guess for the kite's trajectory (states,
    % algebraic variables, and controls) for use in an optimal control problem.
    % It supports generating either a circular or a lemniscate (figure-eight)
    % path in the World frame, then transforms and processes it to derive
    % initial guesses for the system's states, algebraic variables, and controls.
    %
    % Inputs:
    %   symbolic_funcs  - Struct, containing CasADi symbolic functions for:
    %                     .tether_dynamics      : Function describing tether dynamics.
    %                     .tether_alg_equation  : Function for tether algebraic equations.
    %                     .kite_dynamics        : Function describing kite dynamics.
    %                     .winch_dynamics       : Function describing winch dynamics.
    %   params_num      - Struct, numerical parameters of the system.
    %   circRadius      - Scalar, radius of the circular path or characteristic length for lemniscate (m).
    %
    % Outputs:
    %   xi              - Matrix, Initial guess for state variables (N_states x N_points).
    %   zi              - Matrix, Initial guess for algebraic variables (N_alg_vars x N_points).
    %   ui              - Matrix, Initial guess for control variables (N_controls x N_points-1).

    import casadi.*

    % Extract wind direction from parameters
    windDir = params_num.env.wind.direction;
    
    % Extract symbolic functions
    tether_dynamics = symbolic_funcs.tether_dynamics;
    tether_alg_equation = symbolic_funcs.tether_alg_equation;
    kite_dynamics = symbolic_funcs.kite_dynamics;
    winch_dynamics = symbolic_funcs.winch_dynamics;

    % Define the center of the initial trajectory in a plane
    % The trajectory will be generated in a plane that is parallel to the yz-plane
    % of a local frame, then rotated and translated.
    circOffset = 350; % Z-offset from ground (m)
    x_offset = 500;   % X-offset from origin (m)
    el0 = 0 * pi/180; % Initial elevation angle for rotation
    circOrigin_W = [x_offset; 0; circOffset]; % Origin of the circular/lemniscate path in World frame

    % Time discretization parameters
    Tf = params_num.sim.Tf;       % Final time
    N = params_num.sim.Nt;         % Number of time intervals
    h = Tf / N;                    % Time step
    winding_number = params_num.sim.winding_number; % Number of complete loops
    omega = 2 * pi * winding_number / Tf; % Angular velocity for the path

    % Generate time vector for the trajectory points
    t = (0:1/N:1) .* Tf; % Time points from 0 to Tf

    % Initialize position and velocity matrices in a vertical plane
    posKiteVert_W = []; % Position of kite in a vertical plane (World frame)
    velKiteVert_W = []; % Velocity of kite in a vertical plane (World frame)

    % Choose path type (lemniscate or circle)
    if params_num.sim.path_type == 8 % Lemniscate (figure-eight) path
        [POS_W, VEL_W] = build_leminscate(); % Get CasADi functions for lemniscate
        phi_0 = 2 * pi / 4; % Initial phase angle for the lemniscate

        % Populate position and velocity for each time point
        for i = 1:length(t)
            posKiteVert_W = [posKiteVert_W, full(POS_W(t(i), omega, circRadius, 2 * circRadius, phi_0, circOrigin_W))];
            velKiteVert_W = [velKiteVert_W, full(VEL_W(t(i), omega, circRadius, 2 * circRadius, phi_0, circOrigin_W))];
        end
    else % Default to circular path
        [POS_W, VEL_W] = build_circle(); % Get CasADi functions for circle
        phi_0 = 4 * pi / 4; % Initial phase angle for the circle

        % Populate position and velocity for each time point
        for i = 1:length(t)
            posKiteVert_W = [posKiteVert_W, full(POS_W(t(i), omega, circRadius, phi_0, circOrigin_W))];
            velKiteVert_W = [velKiteVert_W, full(VEL_W(t(i), omega, circRadius, phi_0, circOrigin_W))];
        end
    end

    % Rotate the vertical trajectory around the y-axis for a more realistic path.
    % This `R_y` matrix tilts the plane of the circle/lemniscate.
    R_y = [cos(el0) 0 sin(el0); 0 1 0; -sin(el0) 0 cos(el0)];

    % Apply the rotation to positions and velocities
    posKiteInit_W = R_y * posKiteVert_W;
    velKiteInit_W = R_y * velKiteVert_W;

    % Calculate the angular velocity of the body frame relative to the orientation frame
    % (assumed to be a constant rotation corresponding to the path's angular velocity)
    omega_OB_W = R_y * [omega; 0; 0];

    % Define the rotation matrix from World frame 'W' to Orientation frame 'O' (wind frame).
    % This transforms coordinates such that the x-axis aligns with the wind direction.
    M_OW = [cos(windDir), sin(windDir), 0;
            sin(windDir), -cos(windDir), 0;
            0, 0, -1];
    
    % Transform the angular velocity to the Orientation frame
    kiteOmega_OB_O = M_OW * omega_OB_W;
    
    yawResetCounter = 0; % Counter to handle yaw angle wrap-around for continuity

    % Pre-allocate arrays for efficiency
    kiteDist_tau = zeros(1, length(t));
    kiteLat = zeros(1, length(t));
    kiteLong = zeros(1, length(t));
    posKiteInit_O = zeros(3, length(t));
    velKiteInit_O = zeros(3, length(t));
    kiteEuler = zeros(3, length(t));
    kiteOmega_OB_B = zeros(3, length(t));
    velKiteInit_B = zeros(3, length(t));
    winchRotAngle = zeros(1, length(t));
    winchRotVel = zeros(1, length(t));
    tetherAzimuth = zeros(1, length(t));
    tether_elevation = zeros(1, length(t));

    for i = 1:length(t)
        % Calculate kite distance, latitude, and longitude
        kiteDist_tau(i) = norm(posKiteInit_W(:, i));
        kiteLat(i) = asin(posKiteInit_W(3, i) / kiteDist_tau(i));
        kiteLong(i) = atan2(posKiteInit_W(2, i), posKiteInit_W(1, i));
        
        % Transform kite position and velocity to Orientation frame
        posKiteInit_O(:, i) = M_OW * posKiteInit_W(:, i);
        velKiteInit_O(:, i) = M_OW * velKiteInit_W(:, i);

        % Calculate body-frame (B) basis vectors in Orientation frame (O)
        % ex: unit vector along velocity direction
        % ez: unit vector along negative position vector (pointing towards origin)
        % ey: cross product of ez and ex (completing the orthonormal triad)
        ex = velKiteInit_O(:, i) / norm(velKiteInit_O(:, i));
        ez = -posKiteInit_O(:, i) / norm(posKiteInit_O(:, i));
        ey = cross(ez, ex);
        rotMat_OB = [ex ey ez]; % Rotation matrix from Body to Orientation frame

        % Convert rotation matrix to Euler angles (1-2-3 sequence: roll, pitch, yaw)
        if params_num.sim.path_type == 8 % For lemniscate, direct Euler angles conversion
            [kiteEuler(1, i), kiteEuler(2, i), kiteEuler(3, i)] = dcm_to_euler_123(rotMat_OB);
        else % For circular path, handle yaw continuity
            [kiteEuler(1, i), kiteEuler(2, i), kiteEuler_tmp(1, i)] = dcm_to_euler_123(rotMat_OB);
            if i > 1
                % Check for yaw angle wrap-around (crossing +/- pi)
                if abs(kiteEuler_tmp(1, i - 1) - kiteEuler_tmp(1, i)) >= (pi + 0 * (pi / 180))
                    yawResetCounter = yawResetCounter + 1; % Increment counter if wrap-around occurs
                end
            end
            % Apply correction for continuous yaw angle (for counter-clockwise flight)
            kiteEuler(3, i) = kiteEuler_tmp(1, i) - (2 * yawResetCounter) * pi;
        end

        % Transform angular velocity and linear velocity to Body frame
        kiteOmega_OB_B(:, i) = rotMat_OB' * kiteOmega_OB_O;
        velKiteInit_B(:, i) = rotMat_OB' * velKiteInit_O(:, i);
        
        % Calculate winch rotation angle and velocity
        winchRotAngle(1, i) = (kiteDist_tau(i)) / params_num.winch.radius;
        
        % Spherical coordinates transformation matrix for velocity
        sin_lat = sin(kiteLat(i));
        sin_long = sin(kiteLong(i));
        cos_lat = cos(kiteLat(i));
        cos_long = cos(kiteLong(i));
        
        % Rotation matrix from World frame to 'tether-aligned' frame (tau)
        % This frame might have its z-axis aligned with the tether.
        rotMat_tauW = [-sin_lat * cos_long, -sin_lat * sin_long, cos_lat;
                       -sin_long, cos_long, 0;
                       -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];
        
        % Transform kite velocity to tether-aligned frame
        velKiteInit_tau = rotMat_tauW * velKiteInit_W(:, i);
        % Winch rotation velocity derived from the z-component of velocity in tether-aligned frame
        winchRotVel(1, i) = -velKiteInit_tau(3) / params_num.winch.radius;

        % Calculate tether azimuth and elevation angles (used for algebraic variable guess)
        tetherAzimuth(1, i) = asin(posKiteInit_W(2, i) / norm(posKiteInit_W(:, i)));
        tether_elevation(1, i) = acos(posKiteInit_W(3, i) / sqrt(posKiteInit_W(1, i)^2 + posKiteInit_W(3, i)^2));
    end

    % Assemble the initial guess for the state vector (xi)
    % The order of states should match the definition in the problem.
    xi = [winchRotAngle; winchRotVel; velKiteInit_B; kiteOmega_OB_B; kiteEuler; kiteLong; kiteLat; kiteDist_tau];
    
    zi = []; % Initialize algebraic variables vector
    ui = []; % Initialize control variables vector

    % Convert numerical parameters struct to a column vector for CasADi functions
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

    % Optimization options for fsolve (used to solve algebraic equations)
    options = optimoptions('fsolve', 'Algorithm', 'Levenberg-Marquardt', ...
                           'MaxIterations', 3000, 'MaxFunctionEvaluations', 3000, ...
                           'StepTolerance', 1e-10, 'FunctionTolerance', 1e-8, ...
                           'Display', 'none');
    
    tetherForceGuessAtWinch = 9.00e5; % A rough estimation for tether force at winch (N)
    
    % Initial guess for algebraic variables (tether force, elevation, azimuth)
    % using the first point's values.
    z0 = [tetherForceGuessAtWinch; tether_elevation(1); tetherAzimuth(1)];

    % Initial guess for control variables (winch torque, aerodynamic forces/moments)
    % The structure of u0 depends on how controls are defined.
    % Here: [winch_torque; Fx_aero; Fy_aero; Fz_aero; Mx_aero; My_aero; Mz_aero]
    u0 = [2e5; zeros(2, 1); -2e5; zeros(3, 1)]; % Example initialization
    
    tether_tention = 30; % Placeholder value, likely a state or parameter
    numeric_tether_error = 0; % Accumulator for tether algebraic equation error

    % --- Loop through time steps to find initial guesses for algebraic and control variables ---
    for k = 1:N
        % For the first step (k=1), or if ablation study is active, use a specific z0.
        if k == 1
            % Create a MATLAB-compatible function handle for fsolve to solve tether algebraic equations
            % Note: `xi(:,k) - [tether_tention; zeros(13,1)]` suggests that
            % `tether_tention` assumig tether is under tension
            f_tether_matlab = @(Z) full(tether_alg_equation(xi(:, k) - [tether_tention; zeros(13, 1)], Z, paramsColumnVec));
            
            % Re-initialize z0 for each step if ablation study is active, or use previous solution.
            z0 = [tetherForceGuessAtWinch; tether_elevation(k); tetherAzimuth(k)];
            
            % Solve for algebraic variables (tether force at kite, angles)
            [z0_star, ~, ~, ~] = fsolve(f_tether_matlab, z0, options);
            
            % Evaluate tether dynamics with the solved algebraic variables to get errors and forces
            [alg_error, T_kite, ~] = tether_dynamics(xi(:, k) - [tether_tention; zeros(13, 1)], z0_star, paramsColumnVec);
            
            % Store tether forces (winch and kite)
            TetherForces.winch = z0_star(1);
            TetherForces.kite = full(T_kite);
            
            % Solve for control variables (u_star) using the dynamics equations
            u_star = fsolve(@(ui_val)find_u(kite_dynamics, winch_dynamics, xi(:, k + 1), xi(:, k), ui_val, z0_star, paramsColumnVec, h), u0, options);
        end
        
        % Conditional re-initialization of z0 based on ablation study flag (rigid wing)
        if params_num.sim.ablation_study == 1
            z0 = [tetherForceGuessAtWinch; tether_elevation(k); tetherAzimuth(k)];
        else
            z0 = [z0_star(1); tether_elevation(k); tetherAzimuth(k)]; % Use previous tether force estimate
        end
        
        % Solve for algebraic variables using fsolve (hot start from z0 or cold start)
        f_tether_matlab = @(Z) full(tether_alg_equation(xi(:, k) - [tether_tention; zeros(13, 1)], Z, paramsColumnVec));
        [z0_star, ~, ~, ~] = fsolve(f_tether_matlab, z0, options);
        
        % Evaluate tether dynamics to check error and get forces
        [alg_error, T_kite, ~] = tether_dynamics(xi(:, k) - [tether_tention; zeros(13, 1)], z0_star, paramsColumnVec);
        numeric_tether_error = numeric_tether_error + norm(alg_error); % Accumulate error

        % Store tether forces
        TetherForces.winch = z0_star(1);
        TetherForces.kite = full(T_kite);
        
        u0 = u_star; % Use previous control solution as hot start for current step
        % Solve for control variables for the current time step
        u_star = fsolve(@(U)find_u(kite_dynamics, winch_dynamics, xi(:, k + 1), xi(:, k), U, z0_star, paramsColumnVec, h), u0, options);
        
        % Append the solved controls and algebraic variables to their respective lists
        ui = [ui u_star]; 
        zi = [zi (z0_star)]; 
    end

    % Solve for the final algebraic variables at the last state point (k=N+1)
    if params_num.sim.ablation_study == 1
        z0 = [z0_star(1); tether_elevation(N + 1); tetherAzimuth(N + 1)];
    else
        z0 = z0_star; % Use the last solved z_star as hot start
    end
    
    [z0_star, ~, ~, ~] = fsolve(f_tether_matlab, z0, options);
    [alg_error, T_kite, ~] = tether_dynamics(xi(:, end), z0_star, paramsColumnVec);
    zi = [zi (z0_star)]; 
    
    % Report the average error of algebraic equations
    if full(numeric_tether_error) / N > 10.5
        disp('The initial guess might be not appropriate!')
        msg = sprintf('The average error of solving algebraic equations: %d', full(numeric_tether_error) / N);
        disp(msg)
    else
        disp('The total error of solving algebraic equations:')
        disp(full(numeric_tether_error))
    end
    
    % --- Plotting the initialized states, algebraic variables, and controls ---
    
    %  Figure 100: initialized path (for visualization and debugging)
    figure(100)
    plot3(posKiteInit_W(1, 1), posKiteInit_W(2, 1), posKiteInit_W(3, 1), '*b')
    hold on
    plot3(posKiteInit_W(1, 5), posKiteInit_W(2, 5), posKiteInit_W(3, 5), '*r')
    plot3(posKiteInit_W(1, :), posKiteInit_W(2, :), posKiteInit_W(3, :), '-b')
    grid on
    hold on
    xlim([0, 800])
    ylim([-300, 300])
    zlim([0, 700])
    view(30, 21); % Set the view angle for the figure
    
    % Figure 1: States (xi) over time
    figure(1)
    fig_states = gcf;  % Get current figure handle
    set(fig_states, 'Name', 'fig_states', 'units','normalized','Position', [0.1 0.1 0.8 0.8]);

    subplot(7, 1, 1)
    plot(t.', xi(1, :).');
    ylim([min(xi(1, :).') - 1, max(xi(1, :).') + 1])
    title('Winch Rotation Angle');
    ylabel('Angle (rad)');
    grid on
    
    subplot(7, 1, 2)
    plot(t.', xi(2, :).');
    title('Winch Rotation Velocity');
    ylim([min(xi(2, :).') - 1, max(xi(2, :).') + 1])
    ylabel('Velocity (rad/s)');
    grid on
    
    subplot(7, 1, 3)
    plot(t.', xi(3:5, :).');
    title('Kite Velocity in Body Frame');
    ylabel('Velocity (m/s)');
    legend('Vx', 'Vy', 'Vz', 'Location', 'best');
    grid on
    
    subplot(7, 1, 4)
    plot(t.', (180 / pi) .* xi(6:8, :).');
    title('Kite Angular Velocity in Body Frame');
    ylabel('Omega (deg/s)');
    legend('Wx', 'Wy', 'Wz', 'Location', 'best');
    grid on
    
    subplot(7, 1, 5)
    plot(t.', xi(9:11, :).');
    title('Kite Euler Angles');
    ylabel('Angle (rad)');
    legend('Roll', 'Pitch', 'Yaw', 'Location', 'best');
    grid on
    
    subplot(7, 1, 6)
    plot(t.', xi(12:13, :).');
    title('Kite Longitude and Latitude');
    ylabel('Angle (rad)');
    legend('Longitude', 'Latitude', 'Location', 'best');
    grid on
    
    subplot(7, 1, 7)
    plot(t.', xi(14, :).');
    title('Kite Distance from Winch');
    xlabel('Time (s)');
    ylabel('Distance (m)');
    ylim([min(xi(14, :).') - 1, max(xi(14, :).') + 1])
    grid on

    % Figure 2: Algebraic Variables (zi) over time
    figure(2)
    subplot(3, 1, 1)
    plot(t.', zi(1, :).');
    title('Tether Force at Winch');
    ylabel('Force (N)');
    grid on
    
    subplot(3, 1, 2)
    plot(t.', (180 / pi) * zi(2, :).');
    title('Tether Elevation Angle');
    ylabel('Angle (deg)');
    grid on
    
    subplot(3, 1, 3)
    plot(t.', (180 / pi) * zi(3, :).');
    title('Tether Azimuth Angle');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    grid on
    
    % Figure 10: Control Inputs (ui) over time
    figure(10)
    subplot(3, 1, 1)
    stairs(t(1:end-1).', ui(1, :).'); % stairs for control inputs as they are usually piece-wise constant
    title('Winch Motor Torque');
    ylabel('Torque (Nm)');
    grid on
    
    subplot(3, 1, 2)
    stairs(t(1:end-1).', ui(2:4, :).');
    title('Aerodynamic Forces in Body Frame');
    ylabel('Force (N)');
    legend('Fx_a_B', 'Fy_a_B', 'Fz_a_B', 'Location', 'best');
    grid on
    
    subplot(3, 1, 3)
    stairs(t(1:end-1).', ui(5:7, :).');
    title('Aerodynamic Moments in Body Frame');
    xlabel('Time (s)');
    ylabel('Moment (Nm)');
    legend('Mx_a_B', 'My_a_B', 'Mz_a_B', 'Location', 'best');
    grid on
end

%---
% Helper functions for path_init
%---

function ode_eq = find_u(kite_dynamics, winch_dynamics, xi_next, xi, ui, zi, params_num_array, h)
    % FIND_U Calculates the residual of the ODE for solving control inputs.
    %
    % This function computes the difference between the actual next state and
    % the state predicted by the dynamics model, given the current state and
    % an assumed control input. It's used within fsolve to find control inputs
    % that satisfy the discretized dynamics.
    %
    % Inputs:
    %   kite_dynamics   - CasADi function for kite's ODEs.
    %   winch_dynamics  - CasADi function for winch's ODEs.
    %   xi_next         - Vector, State at the next time step.
    %   xi              - Vector, State at the current time step.
    %   ui              - Vector, Current control input guess (winch torque, F_aero, M_aero).
    %   zi              - Vector, Current algebraic variables (tether force, etc.).
    %   params_num_array- Column vector of numerical parameters.
    %   h               - Scalar, Time step (h = Tf / N).
    %
    % Output:
    %   ode_eq          - Vector, Residual of the discretized ODE. This should be
    %                     close to zero if ui is the correct control input.

    % Extract components from the control input vector `ui`
    tether_force_at_winch = zi(1); % Tether force at winch (first algebraic variable)
    u_winch = ui(1);               % Winch motor torque (first control input)
    F_a_B = ui(2:4);               % Aerodynamic forces in body frame (next 3 control inputs)
    M_a_B = ui(5:7);               % Aerodynamic moments in body frame (last 3 control inputs)
    F_p_B = 0;                     % Placeholder for propulsion force in body frame (assumed zero here)

    % Evaluate the right-hand side (RHS) of the Ordinary Differential Equations (ODEs)
    % This combines the winch dynamics and kite dynamics.
    rhs = full(vertcat(winch_dynamics(xi_next, u_winch, tether_force_at_winch, params_num_array), ...
                       kite_dynamics(xi_next, zi, F_a_B, M_a_B, F_p_B, params_num_array)));
    
    % Calculate the residual: xi_next - xi - h * rhs
    % This is the implicit Euler discretization
    % of the ODE: dx/dt = rhs => x_next - x = h * rhs
    ode_eq = xi_next - xi - h .* rhs;
end

%---

function [phi, theta, psi] = dcm_to_euler_123(R)
    % DCM_TO_EULER_123 Converts a Direction Cosine Matrix (DCM) to Euler angles (Z-Y-X sequence).
    %
    % This function converts a 3x3 rotation matrix (Direction Cosine Matrix)
    % into Euler angles following the Z-Y-X (yaw-pitch-roll) convention, also
    % known as 1-2-3 sequence (rotation about X, then Y, then Z).
    %
    % Inputs:
    %   R       - 3x3 matrix, Direction Cosine Matrix (rotation matrix).
    %
    % Outputs:
    %   phi     - Scalar, Roll angle (rotation about X-axis), in radians.
    %   theta   - Scalar, Pitch angle (rotation about Y-axis), in radians.
    %   psi     - Scalar, Yaw angle (rotation about Z-axis), in radians.
    %
    % Reference: Based on standard aerospace conventions for Euler angle extraction.
    % Caution: This conversion has singularities at pitch = +/- 90 degrees.

    phi = atan2(R(3, 2), R(3, 3));   % Roll (rotation about X)
    theta = -asin(R(3, 1));          % Pitch (rotation about Y)
    psi = atan2(R(2, 1), R(1, 1));   % Yaw (rotation about Z)
end

%---

function [POS_W, VEL_W] = build_leminscate()
    % BUILD_LEMINSCATE Creates CasADi functions for a lemniscate (figure-eight) path.
    %
    % This function defines the parametric equations for a lemniscate (Booth's lemniscate)
    % in a 2D plane, then projects and translates it into a 3D World frame.
    % It returns CasADi Function objects for position and velocity, which can
    % be evaluated numerically later.
    %
    % Outputs:
    %   POS_W   - CasADi Function, Takes (t, omega, a, b, phi_0, circOrigin_W)
    %             and returns the 3D position vector in the World frame.
    %   VEL_W   - CasADi Function, Takes (t, omega, a, b, phi_0, circOrigin_W)
    %             and returns the 3D velocity vector in the World frame.

    import casadi.*

    % Declare symbolic variables
    t = SX.sym('t', 1);             % Time
    a = SX.sym('a', 1);             % Lemniscate parameter (scaling factor)
    b = SX.sym('b', 1);             % Lemniscate parameter (related to shape/aspect ratio)
    omega = SX.sym('omega', 1);     % Angular frequency of traversal
    phi_0 = SX.sym('phi_0', 1);     % Initial phase angle
    circOrigin_W = SX.sym('circOrigin_W', 3); % Origin of the lemniscate plane in World frame

    % Normal vector of the plane containing the lemniscate.
    % This assumes the plane passes through the origin of the World frame and
    % is defined by its normal vector pointing from the World origin to circOrigin_W.
    normal_vec = circOrigin_W ./ norm(circOrigin_W);
    
    % Basis vectors in the plane (orthonormal)
    % v1 is perpendicular to normal_vec and lies in the XY plane (if normal_vec_z is non-zero)
    v1 = [-normal_vec(2); normal_vec(1); 0];
    v1 = v1 / norm(v1); % Normalize v1
    v2 = cross(normal_vec, v1); % v2 is perpendicular to both normal_vec and v1

    % Parametric equations for a lemniscate of Booth in 2D (often used for figure-eight shapes)
    % x_2D = a * cos(theta) / (1 + sin(theta)^2)
    % y_2D = b * cos(theta) * sin(theta) / (1 + sin(theta)^2)
    % where theta = omega*t - phi_0
    cos_term = cos(omega * t - phi_0);
    sin_term = sin(omega * t - phi_0);
    denom = 1 + sin_term^2;
    
    x2D = a * cos_term / denom;
    y2D = b * cos_term * sin_term / denom;
    
    % Map the 2D lemniscate coordinates to 3D World frame coordinates
    % The lemniscate is centered at circOrigin_W and lies in the plane spanned by v1 and v2.
    x = circOrigin_W(1) + x2D * v1(1) + y2D * v2(1);
    y = circOrigin_W(2) + x2D * v1(2) + y2D * v2(2);
    z = circOrigin_W(3) + x2D * v1(3) + y2D * v2(3);

    p_k_w = [x; y; z]; % Position vector in World frame
    
    v_k_w = jacobian(p_k_w, t); % Velocity vector (time derivative of position)
    
    % Create CasADi Function objects
    POS_W = Function('POS_W', {t, omega, a, b, phi_0, circOrigin_W}, {p_k_w});
    VEL_W = Function('VEL_W', {t, omega, a, b, phi_0, circOrigin_W}, {v_k_w});
end

%---

function [POS_W, VEL_W] = build_circle()
    % BUILD_CIRCLE Creates CasADi functions for a circular path.
    %
    % This function defines the parametric equations for a circle in a 2D plane,
    % then projects and translates it into a 3D World frame.
    % It returns CasADi Function objects for position and velocity, which can
    % be evaluated numerically later.
    %
    % Outputs:
    %   POS_W   - CasADi Function, Takes (t, omega, circle_r, phi_0, circOrigin_W)
    %             and returns the 3D position vector in the World frame.
    %   VEL_W   - CasADi Function, Takes (t, omega, circle_r, phi_0, circOrigin_W)
    %             and returns the 3D velocity vector in the World frame.

    import casadi.*
    
    % Declare symbolic variables
    t = SX.sym('t', 1);             % Time
    omega = SX.sym('omega', 1);     % Angular frequency of traversal
    circle_r = SX.sym('circle_r', 1); % Radius of the circle
    phi_0 = SX.sym('phi_0', 1);     % Initial phase angle
    circOrigin_W = SX.sym('circOrigin_W', 3); % Origin of the circular plane in World frame

    % Normal vector of the plane containing the circle.
    % Similar to lemniscate, assumes plane passes through World origin to circOrigin_W.
    normal_vec = circOrigin_W ./ norm(circOrigin_W);
    
    % Basis vectors in the plane (orthonormal)
    v1 = [-normal_vec(2); normal_vec(1); 0];
    v1 = v1 / norm(v1); 
    v2 = cross(normal_vec, v1); 

    % Define the circle in the 3D plane.
    x = circOrigin_W(1) + circle_r * (v2(1) * cos(omega * t + phi_0) - v1(1) * sin(omega * t + phi_0));
    y = circOrigin_W(2) + circle_r * (v2(2) * cos(omega * t + phi_0) - v1(2) * sin(omega * t + phi_0));
    z = circOrigin_W(3) + circle_r * (v2(3) * cos(omega * t + phi_0) - v1(3) * sin(omega * t + phi_0));

    p_k_w = [x; y; z]; % Position vector in World frame
    
    v_k_w = jacobian(p_k_w, t); % Velocity vector (time derivative of position)
    
    % Create CasADi Function objects
    POS_W = Function('POS_W', {t, omega, circle_r, phi_0, circOrigin_W}, {p_k_w});
    VEL_W = Function('VEL_W', {t, omega, circle_r, phi_0, circOrigin_W}, {v_k_w});
end