function [phi, theta, psi] = quat_to_euler(q0, qv)
% QUAT_TO_EULER Converts a unit quaternion to Z-Y-X (yaw-pitch-roll) Euler angles.
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
%
% Description:
%   Used only to build an Euler-angle "shadow copy" of a quaternion-based
%   trajectory for the pre-existing Euler-based plotting/animation functions
%   (see others/*.m), which are otherwise left unchanged. Matches the exact
%   convention of the previous dcm_to_euler_123 helper (path_init.m) so the
%   plotted angles are physically identical to the old Euler-native pipeline.
%
% Inputs:
%   q0 - scalar, quaternion scalar part.
%   qv - 3x1, quaternion vector part [q1;q2;q3].
%
% Outputs:
%   phi   - Roll angle (rotation about X), radians.
%   theta - Pitch angle (rotation about Y), radians.
%   psi   - Yaw angle (rotation about Z), radians.

R = quat_to_dcm(q0, qv);

phi   = atan2(R(3,2), R(3,3));   % Roll (rotation about X)
theta = -asin(R(3,1));           % Pitch (rotation about Y)
psi   = atan2(R(2,1), R(1,1));   % Yaw (rotation about Z)

end
