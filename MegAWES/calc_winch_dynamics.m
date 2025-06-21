function dx = calc_winch_dynamics(x, Ft_ground_W, control_moments, params_sym)
    % CALC_WINCH_DYNAMICS Calculates the dynamics of the winch system.
    %
    % This function models the rotational dynamics of a winch, considering
    % the forces exerted by the tether, internal friction, and external
    % control moments. It outputs the time derivatives of the winch's
    % angular position and angular velocity.
    %
    % Inputs:
    %   x             - Current state vector of the winch (2x1):
    %                   x(1) - Angular position of the winch drum (rad).
    %                   x(2) - Angular velocity of the winch drum (rad/s).
    %   Ft_ground_W   - Scalar, Magnitude of the tether tension force at the
    %                   ground station (in the World frame, along the tether direction).
    %                   (Note: This is assumed to be the component that creates torque).
    %   control_moments - Scalar, External control moment applied to the winch drum (Nm).
    %   params_sym    - struct, Symbolic parameters for the winch.
    %
    % Outputs:
    %   dx            - Time derivative of the winch state vector (2x1):
    %                   dx(1) - Angular velocity (omega_w, which is x(2)).
    %                   dx(2) - Angular acceleration (domega_w).

    % Extract system characteristics from the parameters structure
    
    inertia  = params_sym.winch.inertia;   % Rotational inertia of the winch drum
    radius   = params_sym.winch.radius;    % Radius of the winch drum
    friction = params_sym.winch.friction;  % Viscous friction coefficient

    % Extract current angular velocity from the state vector
    omega_w = x(2);

    % Calculate the angular acceleration of the winch drum (domega_w).
    domega_w = (1 / inertia) * ((Ft_ground_W * radius) - (friction * omega_w) - control_moments);

    % Form the output vector `dx` containing the time derivatives of the states
    dx = vertcat(omega_w, domega_w);
end