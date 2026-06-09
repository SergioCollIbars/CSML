function F = shadow_function_cylindrical(R_planet, r_sun, r_sc)
% Cylindrical shadow function.
%
% F = 1 : sunlight
% F = 0 : full eclipse
%
% Inputs:
%   R_planet : occulting body radius [m]
%   r_sun    : planet -> Sun vector [m]
%   r_sc     : planet -> spacecraft vector [m]
%
% All vectors must be in the same frame and same units.

    % Unit vector from planet to Sun
    e_sun = r_sun / norm(r_sun);

    % Projection of spacecraft position along Sun direction
    s_parallel = dot(r_sc, e_sun);

    % Perpendicular distance from Sun-planet line
    r_perp_vec = r_sc - s_parallel * e_sun;
    d_perp = norm(r_perp_vec);

    % Spacecraft must be behind the planet relative to the Sun
    behind_planet = s_parallel < 0;

    % Inside cylindrical shadow
    inside_shadow = behind_planet && d_perp < R_planet;

    if inside_shadow
        F = 0;
    else
        F = 1;
    end
end