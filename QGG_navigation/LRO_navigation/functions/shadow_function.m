function [F] = shadow_function(R_sun, R_planet, r_sun, r_sc)
% Returns F in [0,1]: 1 = Sun fully visible, 0 = fully occulted
% r_sun: planet->Sun (J2000), r_sc: planet->SC (J2000)

    % --- distances
    d_sc  = norm(r_sc);                 % planet->SC
    d_sun = norm(r_sun - r_sc);         % SC->Sun
    
    % --- guard against impossible geometry
    if d_sc <= R_planet
        % you're inside the body (or numerically on surface): Sun is occulted
        F = 0;
        return
    end
    if d_sun <= R_sun
        % you're inside the Sun (won't happen physically here)
        F = 0;
        return
    end

    % --- apparent angular radii on the sky (radians)
    th_s = asin( min(1, R_sun    / d_sun) );  % Sun
    th_e = asin( min(1, R_planet / d_sc ) );  % occulting body

    % --- angular separation between Sun center and body center as seen from SC
    denom = d_sc * d_sun;
    csep  = -(r_sc' * (r_sun - r_sc)) / denom;   % cos(theta_es)
    csep  = max(-1, min(1, csep));
    th_es = acos(csep);

    % --- cases
    % No overlap
    if th_es >= th_s + th_e
        F = 1;
        return
    end

    % Full occultation of the Sun by the body
    if th_e >= th_s && th_es <= (th_e - th_s)
        F = 0;
        return
    end

    % Body completely inside Sun disk.
    % overlap area is just body disk area -> still Sun mostly visible.
    if th_s > th_e && th_es <= (th_s - th_e)
        overlap = pi * th_e^2;
        F = 1 - overlap / (pi * th_s^2);
        return
    end

    % --- partial overlap: use standard two-circle overlap area
    % clamp acos arguments to [-1,1]
    arg1 = (th_es^2 + th_s^2 - th_e^2) / (2 * th_es * th_s);
    arg2 = (th_es^2 + th_e^2 - th_s^2) / (2 * th_es * th_e);
    arg1 = max(-1, min(1, arg1));
    arg2 = max(-1, min(1, arg2));

    a1 = acos(arg1);
    a2 = acos(arg2);

    % overlap area between two circles radius th_s and th_e separated by th_es
    overlap = th_s^2 * a1 + th_e^2 * a2 ...
            - 0.5 * sqrt( max(0, (-th_es + th_s + th_e) * (th_es + th_s - th_e) * ...
                                   (th_es - th_s + th_e) * (th_es + th_s + th_e)) );

    F = 1 - overlap / (pi * th_s^2);
    F = max(0, min(1, F));  % numerical safety
end

