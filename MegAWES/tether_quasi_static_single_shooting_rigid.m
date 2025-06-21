function [eq, T_kite, pj] = tether_quasi_static_single_shooting_rigid(z, pos_and_velocity, tether_length , params_sym, Ns)
% TETHER_QUASI_STATIC_SINGLE_SHOOTING_RIGID Calculates tether forces assuming a rigid-like tether.
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

% This function implements a simplified quasi-static model of the tether,
% treating it as a series of rigid segments with tension propagation,
% but neglecting drag and distributed mass effects. It iterates from
% the kite towards the ground to find the tether tension profile.
% The goal is to determine the initial conditions (tension and angles at
% the kite) that result in the tether's calculated ground end matching
% the actual ground station position.
%
% Inputs:
%   z               - CasADi SX variable (3x1), Algebraic variables representing:
%                     z(1) - Tether tension magnitude at the kite end.
%                     z(2) - Azimuth angle of the tether at the kite end.
%                     z(3) - Elevation angle of the tether at the kite end.
%   pos_and_velocity - CasADi SX variable (6x1), Contains:
%                     pos_and_velocity(1:3) - Kite position vector in the world frame [x; y; z].
%                     pos_and_velocity(4:6) - Kite velocity vector in the world frame [vx; vy; vz].
%   tether_length   - CasADi SX variable (scalar), Total length of the tether.
%   params_sym      - struct, CasADi symbolic parameters.
%   Ns              - scalar (integer), Number of tether segments.
%
% Outputs:
%   eq              - CasADi SX variable (3x1), Algebraic equation representing the
%                     difference between the actual kite position (as the starting point)
%                     and the calculated ground station position. This should be zero.
%   T_kite          - CasADi SX variable (3x1), Tether tension vector at the ground station.
%   pj              - CasADi SX variable (3xNs), Positions of all tether mass points
%                     in the world frame, where pj(:,i) corresponds to the i-th
%                     mass point starting from the ground (index 1) to the kite (index Ns).
%   
%                    This code is an interpretation and adaptation of concepts and calculations
%                    from the 'Tether_QuasiStatic.m' and its sub-function 'Objfunc.m'
%                    from the original MegAWES repository (Copyright 2021 Delft University of
%                    Technology), which is also licensed under the Apache License, Version 2.0.

    import casadi.*
    
    Ls = tether_length/Ns;
    mj = params_sym.tether.rho * Ls;

    E = params_sym.tether.E;
    A = params_sym.tether.A;
    
    p = pos_and_velocity(1:3);
    v = pos_and_velocity(4:6);

    omega_t = cross(p./(norm(p)^2),v);
    
    Tn = z(1);
    theta = z(2);
    phi = z(3);
    
    FT = SX.zeros(3,Ns); % Tension forces
    pj = SX.zeros(3,Ns); % Mass positions
    vj = SX.zeros(3,Ns); % Mass velocities


    % First element from ground station (Ns)
    FT(:,Ns) = Tn.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)];
    pj(:,Ns) = Ls.*[sin(theta)* cos(phi); sin(phi); cos(theta)*cos(phi)];
    vj(:,Ns) = dot(v,p./norm(p)) * (p./norm(p)) + cross(omega_t,pj(:,Ns));
    

    % All other segments and masses except for segment connected to the kite
    for i = Ns:-1:2
        
        FT(:,i-1) = FT(:,i);
 

        l_i_1 = (norm(FT(:,i-1))/(E*A) + 1)*Ls;
        
        pj(:,i-1) = pj(:,i) + l_i_1.*(FT(:,i-1)/norm(FT(:,i-1)));       
        vj(:,i-1) = dot(v,p./norm(p)) * (p./norm(p)) + cross(omega_t,pj(:,i-1));

    end   
    T0 =  FT(:,1);
    l_i_1 = (norm(T0)/(E*A) + 1)*Ls;
    p0 = pj(:,1) + l_i_1.*(T0/norm(T0));
    
    eq = p - p0;
    T_kite = -T0;
end