function dx = calc_ap4_dynamic(x, F_M_a_B, F_p_B, vel_k_W, F_t_W, params_sym)
%Aircraft - 6DoF dynamics
%
% Inputs:
%    v_k_B      - Kinetic velocity in body reference frame
%    DE2019     - Aircraft parameters
%    F_a_B      - Aerodynamic forces acting on kite in body reference frame
%    F_p_B      - Propulsion and brake forces acting on kite in body
%                 reference frame
%    ENVMT      - Environment characteristics
%    F_t_W      - Tether force pulling at the kite
%    Torques_B  - Roll, pitch and yaw moments in body reference frame
%    vel_k_W    - Kinetic velocity in wind reference frame
%    omega_OB_B - Roll, pitch and yaw rates
%    q0,qv      - Quaternion scalar/vector parts (attitude)
%    pos_W      - Kite position in the wind frame W (Cartesian, m)
%
% Outputs:
%    vel_K_B_dot    - Rate of change in kinetic velocity in body reference
%                     frame
%    M_OB           - Rotation matrix from body to inertial reference frame
%    omega_OB_B_dot - Roll, pitch and yaw accellerations
%    q_dot          - Rate of change of the attitude quaternion [q1;q2;q3;q0]
%    pos_W_dot      - Rate of change of the wind-frame position (= vel_k_W)
%
%
% Author: Sebastian Rapp, Ph.D., Revised by: Dylan Eijkelhof, M.Sc.
% Delft University of Technology
% email address: d.eijkelhof@tudelft.nl
% edited by: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Last Modified: June 12, 2025

%-------------	BEGIN CODE	--------------

mass = params_sym.kite.mass;
Jx =  params_sym.kite.Jx;
Jy = params_sym.kite.Jy;
Jz = params_sym.kite.Jz;
Jxz = params_sym.kite.Jxz;


% x(1:2) are winch states [rotation_angle, rotation_velocity]
v_k_B       = x(3:5);
omega_OB_B  = x(6:8);
qv          = x(9:11);
% x(12:14) is the kite position in the wind frame W (Cartesian, meters).
% Its propagation is handled at the bottom of this file; the position states
% themselves do not enter the force/moment balance.
q0          = x(15);

F_a_B       = F_M_a_B(:,1);
T_a_B       = F_M_a_B(:,2);

M_OB = quat_to_dcm(q0, qv);

F_g_O = [0;0; mass * params_sym.env.g];
F_t_O = transformFromWtoO( params_sym.env.wind.direction, F_t_W);

%Forces to body fixed frame
F_g_t_B  = M_OB'*(F_g_O + F_t_O); 
Forces_B = F_a_B + F_g_t_B + vertcat(F_p_B, casadi.DM.zeros(2));

p=omega_OB_B(1);
q=omega_OB_B(2);
r=omega_OB_B(3);



vel_K_B_dot = (1/mass) .* Forces_B - cross(omega_OB_B, v_k_B);

% Quaternion kinematics: q_dot = 0.5 * q (x) [0; omega_OB_B] (Hamilton product).
% Unlike the Euler-angle kinematics this replaces, this has no 1/cos(theta)
% singularity at pitch = +/-90 deg. The unit-norm constraint on [q0;qv] is not
% enforced here (this is just the ODE); it is an NLP equality constraint at
% node 0 (building_NLP.m), carried to every node by the midpoint defects.
q1 = qv(1); q2 = qv(2); q3 = qv(3);

q0_dot = -0.5*(q1*p + q2*q + q3*r);
q1_dot =  0.5*(q0*p + q2*r - q3*q);
q2_dot =  0.5*(q0*q + q3*p - q1*r);
q3_dot =  0.5*(q0*r + q1*q - q2*p);


J = [Jx, 0, Jxz;
     0, Jy, 0;
     Jxz, 0, Jz]; 

omega_OB_B_dot = (J^-1)*(T_a_B - cross(omega_OB_B, mtimes(J,omega_OB_B)));

%% Position propagation
% Cartesian position in the wind frame W: pos_W_dot = vel_k_W. Like the
% quaternion kinematics above, this replaces a singular chart -- the former
% spherical states [long; lat; h_tau] had a 1/(h_tau*cos(lat)) term (blowing
% up at lat = +/-90 deg, exactly the gimbal-lock analogue for position) and
% 1/h_tau divisions. These rows are linear in the states, so the implicit
% midpoint defects in the OCP integrate them without rational-trig error.
pos_W_dot = vel_k_W;

dx = [vel_K_B_dot; omega_OB_B_dot; q1_dot; q2_dot; q3_dot; pos_W_dot; q0_dot];
end


