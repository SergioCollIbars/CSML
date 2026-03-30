function [EARTH_ECI] = ECI2EARTHSUN(r_Earth, r_Sun, r_SC)
    % given the Earth & spacecraft location in the inertial frame. Compute
    % the rotation from the Inertial to Earth pointing frame.
    
    e_vec = r_Sun - r_SC;
    e     = e_vec./vecnorm(e_vec);

    z_B = - (- r_Earth + r_SC)./vecnorm(- r_Earth + r_SC);
    
    y_B = cross(z_B, e)./ vecnorm(cross(z_B, e));
    
    x_B = cross(y_B, z_B);

    EARTH_ECI = [x_B';y_B';z_B'];
end

