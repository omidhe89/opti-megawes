function M_OB = quat_to_dcm(q0, qv)
% QUAT_TO_DCM Builds the body-to-inertial rotation matrix from a unit quaternion.
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
%
% Description:
%   Converts a quaternion [q0; qv] (q0 scalar part, qv = [q1;q2;q3] vector part)
%   into the 3x3 direction cosine matrix M_OB that rotates a vector from the
%   kite body frame 'B' to the ground-fixed inertial frame 'O'. This matches the
%   same body-to-inertial sense as the Z-Y-X Euler DCM previously used in
%   estimations_ap4.m / calc_ap4_dynamic.m. Works transparently for CasADi
%   symbolic types (SX/MX) and plain doubles.
%
% Inputs:
%   q0 - scalar, quaternion scalar part.
%   qv - 3x1, quaternion vector part [q1;q2;q3].
%
% Outputs:
%   M_OB - 3x3, rotation matrix from body frame to inertial frame.

q1 = qv(1); q2 = qv(2); q3 = qv(3);

% Normalizing by the (generally near-1) squared norm keeps M_OB an orthonormal
% rotation matrix even if q drifts slightly off the unit-norm manifold
% (e.g. due to RK4 integration error in the simulation loop).
n2 = q0^2 + q1^2 + q2^2 + q3^2;

M_OB = [ q0^2+q1^2-q2^2-q3^2,   2*(q1*q2-q0*q3),       2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3),       q0^2-q1^2+q2^2-q3^2,   2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2),       2*(q2*q3+q0*q1),       q0^2-q1^2-q2^2+q3^2 ] / n2;

end
