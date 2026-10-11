function q = euler_to_quat(phi, theta, psi)
% EULER_TO_QUAT Converts Z-Y-X (3-2-1, yaw-pitch-roll) Euler angles to a unit quaternion.
%
% Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
%
% Description:
%   Matches the exact convention of dcm_to_euler_123 / quat_to_euler.m:
%   R_OB = Rz(psi)*Ry(theta)*Rx(phi). psi is not wrapped to [-pi,pi] here --
%   it may be a continuously unwrapped yaw spanning multiple winding_number
%   laps -- so the resulting quaternion varies smoothly through the SU(2)
%   double cover (q(psi+2*pi) = -q(psi)) with no branch/sign discontinuity,
%   unlike converting a DCM built from wrapped angles.
%
% Inputs:
%   phi   - Roll angle (rotation about X), radians.
%   theta - Pitch angle (rotation about Y), radians.
%   psi   - Yaw angle (rotation about Z), radians. May be unwrapped/continuous.
%
% Outputs:
%   q - 4x1, quaternion [q0; q1; q2; q3].

qz = [cos(psi/2);   0;              0;              sin(psi/2)];
qy = [cos(theta/2); 0;              sin(theta/2);   0];
qx = [cos(phi/2);   sin(phi/2);     0;              0];

q = quat_hamilton(qz, quat_hamilton(qy, qx));
q = q / norm(q);

end

function q = quat_hamilton(a, b)
    a0 = a(1); av = a(2:4); b0 = b(1); bv = b(2:4);
    q0 = a0*b0 - dot(av, bv);
    qv = a0*bv + b0*av + cross(av, bv);
    q = [q0; qv];
end
