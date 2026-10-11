function q = dcm_to_quat(R)
% DCM_TO_QUAT Converts a body-to-inertial rotation matrix to a unit quaternion.
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
%
% Description:
%   Extracts the quaternion [q0; q1; q2; q3] from a 3x3 rotation matrix using
%   Shepperd's method (branch on the largest of the four squared components),
%   which is numerically robust for any attitude -- no pitch singularity,
%   unlike Euler-angle extraction. The convention is the exact inverse of
%   quat_to_dcm.m: R rotates vectors from the kite body frame 'B' to the
%   inertial frame 'O'.
%
%   The returned quaternion is normalized, so a slightly non-orthonormal
%   input (e.g. the velocity/position triad built in path_init.m) is
%   tolerated silently -- same tolerance as the atan2/asin Euler extraction
%   this replaces. The overall sign of q is arbitrary (+q and -q encode the
%   same rotation); callers that need continuity along a trajectory must
%   enforce it against the previous sample.
%
% Inputs:
%   R - 3x3, rotation matrix (body-to-inertial), near-orthonormal.
%
% Outputs:
%   q - 4x1, unit quaternion [q0; q1; q2; q3] (scalar part first).

tr = R(1,1) + R(2,2) + R(3,3);
[~, idx] = max([tr, R(1,1), R(2,2), R(3,3)]);

switch idx
    case 1 % q0 is the largest component
        s  = 2 * sqrt(1 + tr);
        q0 = 0.25 * s;
        q1 = (R(3,2) - R(2,3)) / s;
        q2 = (R(1,3) - R(3,1)) / s;
        q3 = (R(2,1) - R(1,2)) / s;
    case 2 % q1 is the largest component
        s  = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        q0 = (R(3,2) - R(2,3)) / s;
        q1 = 0.25 * s;
        q2 = (R(1,2) + R(2,1)) / s;
        q3 = (R(1,3) + R(3,1)) / s;
    case 3 % q2 is the largest component
        s  = 2 * sqrt(1 - R(1,1) + R(2,2) - R(3,3));
        q0 = (R(1,3) - R(3,1)) / s;
        q1 = (R(1,2) + R(2,1)) / s;
        q2 = 0.25 * s;
        q3 = (R(2,3) + R(3,2)) / s;
    otherwise % q3 is the largest component
        s  = 2 * sqrt(1 - R(1,1) - R(2,2) + R(3,3));
        q0 = (R(2,1) - R(1,2)) / s;
        q1 = (R(1,3) + R(3,1)) / s;
        q2 = (R(2,3) + R(3,2)) / s;
        q3 = 0.25 * s;
end

q = [q0; q1; q2; q3];
q = q / norm(q);

end