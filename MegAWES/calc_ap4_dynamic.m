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
%    longlat    - Longtitude, latitude
%    EULER      - Euler angles
%    h_tau      - Height in tau reference frame
%
% Outputs:
%    vel_K_B_dot    - Rate of change in kinetic velocity in body reference
%                     frame
%    M_OB           - Rotation matrix from body to inertial reference frame
%    omega_OB_B_dot - Roll, pitch and yaw accellerations
%    longlat_dot    - Rate of change of longtitude, latitude
%    eta_EULER_dot	- Rate of change of euler angles
%    h_tau_dot      - Rate of change of height in tau reference frame
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
EULER       = x(9:11);
longlat     = x(12:13);
h_tau       = x(14); 

F_a_B       = F_M_a_B(:,1);
T_a_B       = F_M_a_B(:,2);


c_psi = cos(EULER(3));
s_psi = sin(EULER(3));

c_theta = cos(EULER(2));
s_theta = sin(EULER(2));

c_phi = cos(EULER(1));
s_phi = sin(EULER(1)); 

M_OB = [ c_psi*c_theta, c_psi*s_theta*s_phi-s_psi*c_phi, c_psi*s_theta*c_phi+s_psi*s_phi;
         s_psi*c_theta, s_psi*s_theta*s_phi+c_psi*c_phi, s_psi*s_theta*c_phi-c_psi*s_phi;
        -s_theta, c_theta*s_phi, c_theta*c_phi];

long = longlat(1);
lat  = longlat(2);

F_g_O = [0;0; mass * params_sym.env.g];
F_t_O = transformFromWtoO( params_sym.env.wind.direction, F_t_W);

%Forces to body fixed frame
F_g_t_B  = M_OB'*(F_g_O + F_t_O); 
Forces_B = F_a_B + F_g_t_B + vertcat(F_p_B, casadi.DM.zeros(2));

p=omega_OB_B(1);
q=omega_OB_B(2);
r=omega_OB_B(3);



vel_K_B_dot = (1/mass) .* Forces_B - cross(omega_OB_B, v_k_B);
c_theta = cos(EULER(2));
s_theta = sin(EULER(2));

c_phi = cos(EULER(1));
s_phi = sin(EULER(1)); 

eta_EULER_dot = [1 s_phi*s_theta/c_theta c_phi*s_theta/c_theta;
                 0 c_phi -s_phi; 
                 0 s_phi/c_theta c_phi/c_theta]*omega_OB_B;


J = [Jx, 0, Jxz;
     0, Jy, 0;
     Jxz, 0, Jz]; 

omega_OB_B_dot = (J^-1)*(T_a_B - cross(omega_OB_B, mtimes(J,omega_OB_B)));
%% Position propagation
sin_lat  = sin(lat);
sin_long = sin(long);
cos_lat  = cos(lat);
cos_long = cos(long);
 
M_tauW = [-sin_lat*cos_long, -sin_lat*sin_long,cos_lat;
          -sin_long, cos_long, 0;
          -cos_lat*cos_long, -cos_lat*sin_long, -sin_lat];

v_k_tau = M_tauW*vel_k_W;

% the abs(h_tau)->h_tau because it's boundry considered to alsways be positive.
long_dot = v_k_tau(2)/( h_tau*cos(lat));  
lat_dot  = v_k_tau(1)/h_tau; 

h_tau_dot = -v_k_tau(3);

longlat_dot = [long_dot;
               lat_dot];
           
dx = [vel_K_B_dot; omega_OB_B_dot; eta_EULER_dot; longlat_dot; h_tau_dot];
end


