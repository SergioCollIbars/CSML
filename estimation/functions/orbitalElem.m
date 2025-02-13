function [alpha] = orbitalElem(r, v, mu)
    % ------------------------------------------------------------------- %
    %                     ORBITAL ELEMENTS FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 19/09/2022

    % Description: Compute the orbital elements at the current time

    % Input:
    %   r: current position vector
    %   v: current velocity vector
    %   mu: gravitational parameter

    % Output:
    %   alpha: orbital elements vector, order:
    %       e: eccentricity vector
    %       h: angular momentum vector
    %       a: semi-major axis
    %       rho: orbital parameter
    %       f: true anomaly
    %       i: orbit inclination angle
    %       Omega: orbit ascending node angle
    %       omega: orbit argument of periapsis
    % ------------------------------------------------------------------- %

    r_u = r / norm(r);        % unitary position vector
    v_u = v / norm(v);        % unitary velocity vector
    
    K = [0, 0, 1]';           % unitary frame vector
    J = [0, 1, 0]';           % unitary frame vector
    I = [1, 0, 0]';           % unitary frame vector

    h = cross(r, v);          % specific angular momentum vector [km^2 / s]
    h_u = h / vecnorm(h);        % unitary h vector.

    n_omega = cross(K, h_u) / vecnorm(cross(K, h_u));  % unitary nodal vector
    n_omegap = cross(h_u, n_omega);

    e = 1 / mu * cross(v, h) - r_u;          % Eccentricity
    e_u = e / vecnorm(e);
    e_up = cross(h_u, e_u);
    
    if(norm(e) < 0)
        warning('Eccentricity value not allowed.')
    end
    
    rho = vecnorm(h)^2 / mu;                 % orbital parameter [km]
    
    a = rho / (1 - norm(e)^2);               % semi-major axis [km]

    i = acos(dot(h_u, K));                   % inclination angle [rad]
    
    c_Omega = dot(n_omega, I);
    s_Omega = dot(n_omega, J);
    
    Omega = atan2(s_Omega, c_Omega);        % Asc node angle [rad]

    c_omega = dot(e, n_omega);
    s_omega = dot(e, n_omegap);

    omega = atan2(s_omega, c_omega);        % Arg of periapsis [rad]

    c_f = dot(r_u, e_u);
    s_f = dot(r_u, e_up);

    f = atan2(s_f, c_f);                    % true anaomaly [rad]

    alpha = [vecnorm(e), vecnorm(h), a, rho, f, i, Omega, omega];

end