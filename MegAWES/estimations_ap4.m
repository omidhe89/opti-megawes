function [pos_O, v_k_W, M_OB] = estimations_ap4(x, params_sym)
% ESTIMATIONS_AP4 Calculates kinematic estimations for the kite in different frames.
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
% 
% Modified by: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 12, 2025
% Last Modified: June 12, 2025
%
% Description:
%   This function computes various kinematic estimations for the kite, including
%   its position and velocity in different reference frames (Ground/Inertial 'O',
%   Wind 'W', Body 'B'). It also calculates the rotation
%   matrix from the body frame to the inertial frame.
%
%   This code is an interpretation and adaptation of concepts and calculations
%   from the 'Estimation' blcok found in the original MegAWES repository
%   (Copyright 2021 Delft University of Technology), which is also licensed
%   under the Apache License, Version 2.0.
%
% Inputs:
%   x          - CasADi SX variable (Nx, 1), State vector of the system:
%                x(1:2)   - Winch states (tether angle, angular velocity)
%                x(3:5)   - Kite body-frame linear velocities [vx; vy; vz]
%                x(6:8)   - Kite body-frame angular velocities [p; q; r]
%                x(9:11)  - Quaternion vector part [q1; q2; q3]
%                x(12:14) - Kite position in the wind frame W (Cartesian, m)
%                x(15)    - Quaternion scalar part q0
%   params_sym - struct, CasADi symbolic parameters
%
% Outputs:
%   pos_O      - CasADi SX variable (3x1), Kite position in the ground-fixed (inertial) frame 'O'.
%   v_k_W      - CasADi SX variable (3x1), Kite kinematic velocity in the wind-fixed frame 'W'.
%   M_OB       - CasADi SX variable (3x3), Rotation matrix from Body frame 'B' to Inertial frame 'O'.
%
%   NOTE: the W-frame Cartesian position is NOT an output -- it is the state
%   x(12:14) itself; index the state vector directly where it is needed.
%
% See also: build_dynamics, system_vars_params, transformFromWtoO

% Extract wind direction from symbolic parameters.
wind_direction = params_sym.env.wind.direction;

% Extract individual state components from the input state vector 'x'.
% These correspond to the AP4 (Aircraft Part 4) model states.
v_k_B       = x(3:5);  % Kite kinetic velocity in body frame (vx, vy, vz)
qv          = x(9:11); % Quaternion vector part [q1; q2; q3]
pos_W       = x(12:14);% Kite position in the wind frame W (Cartesian, m)
q0          = x(15);   % Quaternion scalar part

% --- Compute Rotation Matrix from Body Frame to Inertial Frame (M_OB) ---
% This matrix transforms vectors from the kite's body-fixed coordinate
% system to the ground-fixed inertial coordinate system (NED - North-East-Down).
% Built from the attitude quaternion [q0; qv]; matches the same body-to-inertial
% sense as the Z-Y-X (yaw-pitch-roll) Euler DCM this replaced.
M_OB = quat_to_dcm(q0, qv);

% --- Compute Rotation Matrix from Wind Frame to Ground-Fixed Frame (M_WO) ---
% This matrix rotates vectors from a wind-aligned frame (often assumed to be
% aligned with the base wind direction) to the ground-fixed inertial frame.
% The transformation typically involves a rotation around the vertical axis.
M_WO = [cos(wind_direction), sin(wind_direction), 0;
        sin(wind_direction), -cos(wind_direction), 0; 
        0,                   0,                   -1]; % Note: -1 for Z-axis inversion from wind frame to inertial (NED)

% --- Calculate Kite Position in Ground-Fixed (Inertial) Frame (pos_O) ---
% Transforms the position from the wind frame to the ground-fixed frame.
pos_O = M_WO' * pos_W; % Transpose of M_WO is used for inverse transformation.

% --- Calculate Kite Velocity in Ground-Fixed (Inertial) Frame (v_k_O) ---
% Transforms the kite's velocity from its body frame to the ground-fixed frame.
v_k_O = M_OB * v_k_B;

% --- Calculate Kite Velocity in Wind Frame (v_k_W) ---
% Transforms the kite's velocity from the ground-fixed frame to the wind frame.
v_k_W = M_WO * v_k_O;

end % End of function estimations_ap4