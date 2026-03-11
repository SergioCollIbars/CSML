function [EARTH_ECI] = ECI2EARTH(r_Earth,v_Earth, r_SC, v_SC)
    % given the Earth & spacecraft location in the inertial frame. Compute
    % the rotation from the Inertial to Earth pointing frame.
    
    k = [0;0;1];

    r = r_SC - r_Earth;  
    % v = v_SC - v_Earth;

    x_B = - r./vecnorm(r);
    y_B = cross(k, x_B);
    z_B = cross(x_B, y_B);

    EARTH_ECI = [x_B';y_B';z_B'];
end

