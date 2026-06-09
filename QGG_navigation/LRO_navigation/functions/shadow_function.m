function F = shadow_function(R_sun, R_planet, r_sun, r_sc)
%SHADOW_FUNCTION Fraction of solar disk visible from spacecraft.
%
% F = 1 : Sun fully visible
% F = 0 : Sun fully occulted by planet
%
% Inputs:
%   R_sun    : Sun radius [same units as vectors]
%   R_planet : occulting body radius [same units as vectors]
%   r_sun    : planet -> Sun vector
%   r_sc     : planet -> spacecraft vector
%
% Both vectors must be expressed in the same frame.

    % Distances
    d_sc  = norm(r_sc);             % planet -> spacecraft
    rho_s = r_sun - r_sc;           % spacecraft -> Sun
    d_sun = norm(rho_s);            % spacecraft -> Sun

    % Default: no shadow if geometry is invalid
    if d_sc <= R_planet || d_sun <= R_sun
        F = NaN;
        return
    end

    % Unit vectors as seen from spacecraft
    u_planet = -r_sc / d_sc;        % spacecraft -> planet
    u_sun    = rho_s / d_sun;       % spacecraft -> Sun

    % Apparent angular radii
    theta_p = asin(R_planet / d_sc);
    theta_s = asin(R_sun    / d_sun);

    % Angular separation between apparent centers
    cos_sep = dot(u_planet, u_sun);
    cos_sep = max(-1, min(1, cos_sep));
    theta_ps = acos(cos_sep);

    % No overlap: Sun fully visible
    if theta_ps >= theta_p + theta_s
        F = 1;
        return
    end

    % Planet fully covers Sun: total eclipse
    if theta_p >= theta_s && theta_ps <= theta_p - theta_s
        F = 0;
        return
    end

    % Sun fully contains planet disk: annular/partial occultation
    if theta_s > theta_p && theta_ps <= theta_s - theta_p
        occulted_area = pi * theta_p^2;
        F = 1 - occulted_area / (pi * theta_s^2);
        return
    end

    % Partial overlap between two apparent disks
    r1 = theta_s;    % apparent Sun radius
    r2 = theta_p;    % apparent planet radius
    d  = theta_ps;   % angular separation

    arg1 = (d^2 + r1^2 - r2^2) / (2*d*r1);
    arg2 = (d^2 + r2^2 - r1^2) / (2*d*r2);

    arg1 = max(-1, min(1, arg1));
    arg2 = max(-1, min(1, arg2));

    area_overlap = ...
        r1^2 * acos(arg1) + ...
        r2^2 * acos(arg2) - ...
        0.5 * sqrt(max(0, ...
        (-d + r1 + r2) * ...
        ( d + r1 - r2) * ...
        ( d - r1 + r2) * ...
        ( d + r1 + r2)));

    F = 1 - area_overlap / (pi * theta_s^2);

    % Numerical cleanup
    F = max(0, min(1, F));
end

