function [SUN_ECI] = ECI2SUN(r_Sun,r_SC)
    % given the Sun & spacecraft location in the inertial frame. Compute
    % the rotation from the Inertial to Sun pointing frame.
    
    k = [0;0;1];
    x_B = (r_Sun - r_SC)./vecnorm(r_Sun - r_SC);

    y_B = cross(k, x_B)./vecnorm(cross(k, x_B));

    z_B = cross(x_B, y_B);
    
    SUN_ECI = [x_B';y_B';z_B'];
end

