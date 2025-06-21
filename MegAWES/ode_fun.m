function dx = ode_fun(x, u, vel_k_W, TetherForces, Params)
    u_w = u(1);
    u_ap4 = u(2:8);
    dx_w = calc_winch_dynamics(x(1:2), TetherForces.winch, u_w, Params);
    dx_ap4 = calc_ap4_dynamic(x(3:14), u_ap4, vel_k_W, TetherForces.kite, Params);
    dx = vertcat(dx_w, dx_ap4);

end