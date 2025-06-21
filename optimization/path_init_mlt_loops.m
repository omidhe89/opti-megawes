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
function [xi, zi, ui] = path_init_mlt_loops(symbolic_funcs, params_sym, params_num, circRadius)

    import casadi.*
    windDir = params_num.env.wind.direction;
    baseWindSpeed = params_num.env.wind.base_velocity;

    tether_dynamics = symbolic_funcs.tether_dynamics;
    tether_alg_equation = symbolic_funcs.tether_alg_equation;
    kite_dynamics = symbolic_funcs.kite_dynamics;
    winch_dynamics = symbolic_funcs.winch_dynamics;
    
    Tf = params_num.sim.Tf;
    N = params_num.sim.Nt;
    h = Tf/N;
    winding_number = params_num.sim.winding_number;
    points_per_loop = N / winding_number;

    % Define time for traction phase
    T_traction = ((winding_number - 0.5) / winding_number) * Tf;
    N_traction = points_per_loop * (winding_number - 0.5);
%     t_trac = linspace(0, T_traction, points_per_loop * (winding_number -1));
    t_trac = (0:1/N_traction:1) * T_traction;
    omega = (2 * pi * (winding_number - 0.5)) / T_traction;

    posKiteVert_W = [];
    velKiteVert_W = [];
    if params_num.sim.path_type == 8
        [POS_W, VEL_W] = build_tract_leminscate();
    else
        [POS_W, VEL_W] = build_tract_circle();
    end
    phi_0 = 0*pi/4;

    % lemniscate of booth
    circOffset = [450; 0; 100];
%     circRadius = 100;
     
     for i=1:length(t_trac)
        posKiteVert_W = [posKiteVert_W, full(POS_W(t_trac(i), omega, circRadius, phi_0, circOffset))];
        velKiteVert_W = [velKiteVert_W, full(VEL_W(t_trac(i), omega, circRadius, phi_0, circOffset))];
     end
     t_retrac = (0:1/points_per_loop:1) * (Tf - T_traction) + T_traction;
%      t_retrac = linspace(T_traction, Tf, points_per_loop);
     [POS_W, VEL_W] = build_retract(T_traction, Tf- T_traction, posKiteVert_W(:,1), velKiteVert_W(:,1), posKiteVert_W(:,end), velKiteVert_W(:,end));
     for i=1:length(t_retrac)
        posKiteVert_W = [posKiteVert_W, full(POS_W(t_retrac(i)))];
        velKiteVert_W = [velKiteVert_W, full(VEL_W(t_retrac(i)))];
     end
     % rotate vertical trajcetory around y-axis for more realistic initial path
    el0 = -30 * pi/180; 
    R_y = [cos(el0) 0 sin(el0);  0 1 0; -sin(el0) 0 cos(el0)];

    posKiteInit_W = R_y * posKiteVert_W;
    velKiteInit_W = R_y * velKiteVert_W;

    omega_OB_W = R_y *[omega; 0; 0];
    
    M_OW = [cos(windDir), sin(windDir), 0;
            sin(windDir), -cos(windDir), 0;
            0, 0, -1];
    t = (0:1/N:1)*Tf;
    kiteOmega_OB_O = M_OW * omega_OB_W;
    yawResetCounter = 0; % to make yaw rotation continuous
    for i=1:length(t)
         
        kiteDist_tau(i) = norm(posKiteInit_W(:,i),2);
        kiteLat(i) = asin(posKiteInit_W(3,i) / kiteDist_tau(1,i));
        kiteLong(i) = atan2(posKiteInit_W(2,i),posKiteInit_W(1,i));    
        
        posKiteInit_O(:,i) = M_OW * posKiteInit_W(:,i);
        velKiteInit_O(:,i) = M_OW * velKiteInit_W(:,i);


        ex = velKiteInit_O(:,i)/norm(velKiteInit_O(:,i));
        ez = -posKiteInit_O(:,i)/ norm(posKiteInit_O(:,i));
        ey = cross(ez,ex);
        rotMat_OB = [ex ey ez];


        if params_num.sim.path_type == 8
            [kiteEuler(1,i), kiteEuler(2,i), kiteEuler(3,i)] = dcm_to_euler_123(rotMat_OB);
        
        else
           
            [kiteEuler(1,i), kiteEuler(2,i), kiteEuler_tmp(1,i)] = dcm_to_euler_123(rotMat_OB);
            if i > 1 
                if abs(kiteEuler_tmp(1,i-1) - kiteEuler_tmp(1,i)) >= (pi - 10*(pi/180))
                     yawResetCounter = yawResetCounter + 1;
                else
                    kiteEuler(3,i) = kiteEuler_tmp(1,i);
                end

            else
                kiteEuler(3,i) = kiteEuler_tmp(1,i);
            end
%               kiteEuler(3,i)  = kiteEuler_tmp(1,i) +  ( 2 * yawResetCounter)*pi; % fly CW
                kiteEuler(3,i)  = kiteEuler_tmp(1,i) -  (2*yawResetCounter)*pi;    % fly CCW

        end

        kiteOmega_OB_B(:,i) = rotMat_OB' * kiteOmega_OB_O;
        velKiteInit_B(:,i) = rotMat_OB' * velKiteInit_O(:,i);
        
        winchRotAngle(1,i) = (kiteDist_tau(i))/params_num.winch.radius;
        sin_lat  = sin(kiteLat(i));
        sin_long = sin(kiteLong(i));
        cos_lat  = cos(kiteLat(i));
        cos_long = cos(kiteLong(i));
 
        rotMat_tauW = [-sin_lat*cos_long, -sin_lat*sin_long,cos_lat;
                -sin_long, cos_long, 0;
                -cos_lat*cos_long, -cos_lat*sin_long, -sin_lat];

        velKiteInit_tau = rotMat_tauW * velKiteInit_W(:,i);
        winchRotVel(1,i) = -velKiteInit_tau(3) / params_num.winch.radius;


        % azithmuth and elevation angles to start finding algebric variables
        tetherAzimuth(1,i) = asin(posKiteInit_W(2,i)./norm(posKiteInit_W(:,i)));
        tether_elevation(1,i) = acos(posKiteInit_W(3,i)./sqrt(posKiteInit_W(1,i).^2 + posKiteInit_W(3,i).^2));

    end

    xi = [winchRotAngle; winchRotVel; velKiteInit_B; kiteOmega_OB_B; kiteEuler; kiteLong; kiteLat; kiteDist_tau];     % - 25 * ones(1, N+1)

    figure(100)
    plot3(posKiteInit_W(1,1),posKiteInit_W(2,1),posKiteInit_W(3,1),'*r')
    hold on
    hold on
    plot3(posKiteInit_W(1,5),posKiteInit_W(2,5),posKiteInit_W(3,5),'*r')
    plot3(posKiteInit_W(1,:),posKiteInit_W(2,:),posKiteInit_W(3,:),'-b')
    grid on
    hold on
    xlim([0, 800])
    ylim([-300, 300])
    zlim([0, 700])
    view(30,21); % set the view angle for figure
    
    zi = [];
    ui = [];

    tetherForceGuessAtWinch = 1.5e6; % a rough estimation about tether force value! 

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
    options = optimoptions('fsolve','Algorithm', 'Levenberg-Marquardt', 'MaxIterations',3000, 'MaxFunctionEvaluations',3000,'StepTolerance',1e-10,'FunctionTolerance',1e-8,'Display','none');
    

%     zi = [zi z0];  
    u0 = [10e5; zeros(2,1); -5e5 ;zeros(3,1)];
    tether_tention = 0;
    numeric_tether_error = 0;
%     ui = [ui u0];
    for k = 1:N
        
        if k == 1
            % Create a MATLAB-compatible function handle for fsolve
            f_tether_matlab = @(zi) full(tether_alg_equation(xi(:,k), zi, paramsColumnVec));
            z0 = [tetherForceGuessAtWinch; tether_elevation(k); tetherAzimuth(k)];  
            [z0_star, fval, exitflag, output] = fsolve(f_tether_matlab, z0, options);
            [alg_error, T_kite, tether_pos] = tether_dynamics(xi(:,k), z0_star, paramsColumnVec);
            TetherForces.winch = z0_star(1);
            TetherForces.kite = T_kite;
%             u0 = zeros(7,1);
            u_star = fsolve(@(ui)find_u(kite_dynamics, winch_dynamics, xi(:,k+1), xi(:,k), ui, z0_star, paramsColumnVec, h), u0, options);
            zi = [zi (z0_star)];         
%             ui = [ui u_star];
        end
        if params_num.sim.ablation_study == 1
            z0_star = [zi(1,end); tether_elevation(k+1); tetherAzimuth(k+1)];
        else
%           Create a MATLAB-compatible function handle for fsolve
            f_tether_matlab = @(Z) full(tether_alg_equation(xi(:, k+1), Z, paramsColumnVec));
            [z0_star, fval, exitflag, output] = fsolve(f_tether_matlab, [tetherForceGuessAtWinch; tether_elevation(k+1); tetherAzimuth(k+1)], options);
            [alg_error, T_kite, tether_pos] = tether_dynamics(xi(:, k+1), z0_star, paramsColumnVec);
            numeric_tether_error = numeric_tether_error + norm(alg_error); % Accumulate error
        end
        TetherForces.winch = z0_star(1);
        TetherForces.kite = T_kite;
        u0 = u_star;
        u_star = fsolve(@(U)find_u(kite_dynamics, winch_dynamics, xi(:,k+1), xi(:,k), U, z0_star, paramsColumnVec, h), u0, options);
        ui = [ui u_star];
        zi = [zi (z0_star)];
    end
    
    % Report the average error of algebraic equations
    if full(numeric_tether_error) / N > 10.5
        disp('The initial guess might be not appropriate!')
        msg = sprintf('The average error of solving algebraic equations: %d', full(numeric_tether_error) / N);
        disp(msg)
    else
        disp('The total error of solving algebraic equations:')
        disp(full(numeric_tether_error)/ N)
    end
    
    
    figure(1)
    subplot(7,1,1)
    plot(t.', xi(1,:).');
    hold on
    grid on
    subplot(7,1,2)
    plot(t.', xi(2,:).');
    hold on
    grid on
    subplot(7,1,3)
    plot(t.', xi(3:5,:).');
    hold on
    grid on
    subplot(7,1,4)
    plot(t.', (180/pi) .* xi(6:8,:).');
    hold on
    grid on
    subplot(7,1,5)
    plot(t.', (180/pi) .* xi(9:11,:).');
    hold on
    grid on
    subplot(7,1,6)
    plot(t.', xi(12:13,:).');
    hold on
    grid on
    subplot(7,1,7)
    plot(t.', xi(14,:).');
    hold on
    grid on

    
    figure(2)
    subplot(3,1,1)
    plot(t.', zi(1,:).');
    hold on
    grid on
    subplot(3,1,2)
    plot(t.', (180/pi) * zi(2,:).');
    hold on
    grid on
    subplot(3,1,3)
    plot(t.', (180/pi) * zi(3,:).');
    hold on
    grid on
   
    

    figure(10)
    subplot(3,1,1)
    stairs(t(1:end-1).', ui(1,:).');
%     plot(t.', zeros(1, length(t)));
    ylabel('M_winch')
    subplot(3,1,2)
    stairs(t(1:end-1).', ui(2:4,:).');
    ylabel('F_a_B')
    subplot(3,1,3)
    stairs(t(1:end-1).', ui(5:7,:).');
    ylabel('M_a_B')
end

function [POS_W, VEL_W] = build_tract_leminscate()
    import casadi.*
    t = SX.sym('t',1);
    radius = SX.sym('r',1);
    a = 2*radius;
    b = radius;
    omega = SX.sym('omega',1);
    phi_0 = SX.sym('phi_0',1);
    offset =  SX.sym('offset',3);
    % Lemniscate trajectory (shifted to start cleanly)
    x =  offset(1) + t / omega;
    y =  offset(2) + a * cos(omega * t + phi_0);
    z =  offset(3) + b * sin(omega * t + phi_0) .* cos(omega * t + phi_0);
    p_k_w = [x;y;z];
    v_k_w = jacobian(p_k_w, t);
    POS_W = Function('POS_W', {t, omega, radius, phi_0, offset}, {p_k_w});
    VEL_W = Function('VEL_W', {t, omega, radius, phi_0, offset}, {v_k_w});
end  
function [POS_W, VEL_W] = build_tract_circle()
    import casadi.*
    t = SX.sym('t',1);
    radius = SX.sym('r',1);
    omega = SX.sym('omega',1);
    phi_0 = SX.sym('phi_0',1);
    offset =  SX.sym('offset',3);
    % Lemniscate trajectory (shifted to start cleanly)
    x = offset(1) + 2 * t / omega;
    y = offset(2) + radius * cos(omega * t + phi_0);
    z = offset(3) + radius * sin(omega * t + phi_0); % Swapped
    p_k_w = [x;y;z];
    v_k_w = jacobian(p_k_w, t);
    POS_W = Function('POS_W', {t, omega, radius, phi_0, offset}, {p_k_w});
    VEL_W = Function('VEL_W', {t, omega, radius, phi_0, offset}, {v_k_w});
end 
function [POS_W, VEL_W] = build_retract(T_trac, T_retrac, p_0_W, v_0_W, p_end_W, v_end_W)
    import casadi.*
    t = SX.sym('t',1);
    s = (t-T_trac)/T_retrac;

    % Hermite basis functions for smooth transition
    h00 = (1 + 2*s) .* (1 - s).^2;
    h10 = s .* (1 - s).^2;
    h01 = s.^2 .* (3 - 2*s);
    h11 = s.^2 .* (s - 1);


    % Interpolated return trajectory
    x_return = h00 * p_end_W(1) + h10 * (v_end_W(1) * T_retrac) + h01 * p_0_W(1) + h11 * (v_0_W(1) * T_retrac);
    y_return = h00 * p_end_W(2) + h10 * (v_end_W(2) * T_retrac) + h01 * p_0_W(2) + h11 * (v_0_W(2) * T_retrac);
    z_return = h00 * p_end_W(3) + h10 * (v_end_W(3) * T_retrac) + h01 * p_0_W(3) + h11 * (v_0_W(3) * T_retrac);
    p_k_w = [x_return; y_return; z_return];
    v_k_w = jacobian(p_k_w, t);
    POS_W = Function('POS_W', {t}, {p_k_w});
    VEL_W = Function('VEL_W', {t}, {v_k_w});
end
function ode_eq = find_u(kite_dynamics, winch_dynamics, xi_next, xi, ui, zi, params_num_array, h)
    
    tether_force_at_winch = zi(1);
    u_winch = ui(1);
    F_a_B = ui(2:4);
    M_a_B = ui(5:7);

    F_p_B = 0;
    rhs = full(vertcat(winch_dynamics(xi_next, u_winch, tether_force_at_winch, params_num_array), kite_dynamics(xi_next, zi, F_a_B, M_a_B, F_p_B, params_num_array)));
    ode_eq = xi_next - xi - h.*rhs;
end
function [phi, theta, psi] = dcm_to_euler_123(R)
 
    phi = atan2(R(3, 2), R(3, 3));
    theta = -asin(R(3, 1));
    psi = atan2(R(2, 1), R(1, 1));
end