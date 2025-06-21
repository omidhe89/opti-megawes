% TETHER_QUASI_STATIC_SINGLE_SHOOTING Solves for tether shape and forces using a quasi-static model.
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
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
%
function [eq, T_kite, pj] = tether_quasi_static_single_shooting(z, pos_and_velocity, tether_length, params_sym, Ns)

% Description:
%   This function implements a quasi-static model of the tether, dividing it
%   into 'Ns' segments. It calculates the tension and position of each
%   segment's mass point, starting from the kite end and iterating towards
%   the ground station. The goal is to find the initial conditions (tension
%   and angles at the ground) that result in the tether end matching the
%   known kite position and velocity. This is formulated as an algebraic
%   equation `eq = p - p0`, where `p` is the actual kite position and `p0`
%   is the calculated tether end position from the ground.
%
%   This code is an interpretation and adaptation of concepts from the
%   'Tether_QuasiStatic.m' and its sub-function 'Objfunc.m' from the
%   original MegAWES repository (Copyright 2021 Delft University of Technology),
%   which is also licensed under the Apache License, Version 2.0.
%
% Inputs:
%   z               - CasADi SX variable (3x1), Algebraic variables:
%                     z(1) - Tether tension at the ground station (initial guess).
%                     z(2) - Azimuth angle of the tether at the ground station (initial guess).
%                     z(3) - Elevation angle of the tether at the ground station (initial guess).
%   pos_and_velocity - CasADi SX variable (6x1), Contains:
%                     pos_and_velocity(1:3) - Kite position vector in the world frame [x; y; z].
%                     pos_and_velocity(4:6) - Kite velocity vector in the world frame [vx; vy; vz].
%   tether_length   - CasADi SX variable (scalar), Total length of the tether.
%   params_sym      - struct, CasADi symbolic parameters.
%   Ns              - scalar (integer), Number of tether segments.
%
% Outputs:
%   eq              - CasADi SX variable (3x1), Algebraic equation representing the
%                     difference between the actual kite position and the
%                     calculated tether end position from the ground. This
%                     should be zero at equilibrium.
%   T_kite          - CasADi SX variable (3x1), Tether tension vector at the kite (segment Ns).
%   pj              - CasADi SX variable (3xNs), Positions of all tether mass points
%                     in the world frame, where pj(:,i) corresponds to the i-th
%                     mass starting from the ground (index 1) to the kite (index Ns).
%
% See also: build_dynamics, log_wind  
    
    import casadi.*
    base_wind_speed = params_sym.env.wind.base_velocity;
    Ls = tether_length/Ns;
    mj = params_sym.tether.rho * Ls;

    E = params_sym.tether.E;
    A = params_sym.tether.A;
    g = params_sym.env.g;
    
    p = pos_and_velocity(1:3);
    v = pos_and_velocity(4:6);

    omega_t = cross(p./(norm(p)^2),v);
    
    Tn = z(1);
    theta = z(2);
    phi = z(3);
    
    FT = SX.zeros(3,Ns); % Tension forces
    pj = SX.zeros(3,Ns); % Mass positions
    vj = SX.zeros(3,Ns); % Mass velocities
    aj = SX.zeros(3,Ns); % Mass accellerations
    windVel = SX.zeros(3,Ns); % mass wind speed
    Fd = SX.zeros(3,Ns); % drag forces

    % First element from ground station (Ns)
    FT(:,Ns) = Tn.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)];
    pj(:,Ns) = Ls.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)];
    vj(:,Ns) = dot(v,p./norm(p)) .* (p./norm(p)) + cross(omega_t,pj(:,Ns));
    aj(:,Ns) = cross(omega_t,cross(omega_t,pj(:,Ns)));
    
    % Drag calculation first element
%     H_rho = 8550; % in meters (8.55 km)
%     rho_factor = exp(pj(3,Ns) / (2 * H_rho));
    rho_factor  = 1;
    windVel(:,Ns) = [0; 0; log_wind(pj(3,Ns), base_wind_speed, params_sym)];
    v_a_p = vj(:,Ns) - windVel(:,Ns);
    v_a_p_t = dot((pj(:,Ns))/norm(pj(:,Ns)),v_a_p)*((pj(:,Ns))/norm(pj(:,Ns)));
    v_a_p_n = v_a_p - v_a_p_t;
    rho = params_sym.env.wind.rho * rho_factor;
    Fd(:,Ns) = -0.5 * rho* Ls * params_sym.tether.d * params_sym.tether.cd * ...
            norm(v_a_p_n) * v_a_p_n; % particle drag


    % All other segments and masses except for segment connected to the kite
    for i = Ns:-1:2
        if i == Ns
            FT(:,i-1) = (mj+0.5*mj).*aj(:,i) + FT(:,i) - Fd(:,i) + ...
                [0;0;(mj+0.5*mj)*g];
        else
            FT(:,i-1) = mj.*aj(:,i) + FT(:,i) - Fd(:,i) + [0;0;mj*g];
        end

        l_i_1 = (norm(FT(:,i-1))/(E*A) + 1)*Ls;
        
        pj(:,i-1) = pj(:,i) + l_i_1.*(FT(:,i-1)/norm(FT(:,i-1)));       
        vj(:,i-1) = dot(v,p./norm(p)) .* (p./norm(p)) + cross(omega_t,pj(:,i-1));
        aj(:,i-1) = cross(omega_t,cross(omega_t,pj(:,i-1)));
        
        % Drag calculation
%         rho_factor = exp(pj(3,i) / (2 * H_rho));
        rho_factor = 1;
        windVel(:,i) = [0; 0; log_wind(pj(3,i), base_wind_speed, params_sym)];
        v_a_p = vj(:,i) - windVel(:,i);
        v_a_p_t = dot((pj(:,i-1)-pj(:,i))/norm(pj(:,i-1)-pj(:,i)),v_a_p) * ...
                ((pj(:,i-1)-pj(:,i))/norm(pj(:,i-1)-pj(:,i)));
        v_a_p_n = v_a_p - v_a_p_t;
        rho = params_sym.env.wind.rho * rho_factor;
        Fd(:,i-1) = -0.5 * rho * Ls * params_sym.tether.d * params_sym.tether.cd * ...
                norm(v_a_p_n) * v_a_p_n; % particle drag
    end   
    T0 = (mj+0.5*mj).*aj(:,1) + FT(:,1) - Fd(:,1) + [0;0;(mj+0.5*mj)*g];
    l_i_1 = (norm(T0)/(E*A) + 1)*Ls;
    p0 = pj(:,1) + l_i_1.*(T0/norm(T0));
    
    eq = p - p0;
    T_kite = -T0;
end