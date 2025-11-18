function rates = j2_secular_rates(a, e, inc, mu, Re, J2)
    % J2 secular rates for RAAN (Omega), argument of perigee (omega), and mean anomaly (M).
    % Inputs:
    %   a   - semi-major axis [m]
    %   e   - eccentricity [-]
    %   inc - inclination [rad]
    %   mu  - gravitational parameter [m^3/s^2]
    %   Re  - equatorial radius [m]
    %   J2  - second zonal harmonic [-]
    %
    % Vectorized: a, e, inc can be scalars or same-sized arrays.
    %
    % Outputs (struct 'rates'):
    %   .Omega_rad_s  - dOmega/dt [rad/s]
    %   .omega_rad_s  - domega/dt [rad/s]
    %   .M_rad_s      - dM/dt [rad/s]
    %   .Omega_deg_day, .omega_deg_day, .M_deg_day - same, in [deg/day]

    arguments
        a   {mustBeReal, mustBePositive}
        e   {mustBeReal}
        inc {mustBeReal}
        mu  (1,1) {mustBeReal, mustBePositive}
        Re  (1,1) {mustBeReal, mustBePositive}
        J2  (1,1) {mustBeReal}
    end

    % Basic quantities
    n  = sqrt(mu ./ (a.^3));                   % mean motion [rad/s]
    p  = a .* (1 - e.^2);                      % semilatus rectum [m]
    cI = cos(inc);
    srt = sqrt(max(0, 1 - e.^2));              % guard tiny negatives from rounding

    % Secular J2 rates (Lagrange planetary equations, orbit-averaged)
    dOmega = -(3/2) .* n .* J2 .* (Re./p).^2 .* cI;                       % [rad/s]
    domega =  (3/4) .* n .* J2 .* (Re./p).^2 .* (5*cI.^2 - 1);            % [rad/s]
    dM     =  n + (3/4) .* n .* J2 .* (Re./p).^2 .* srt .* (3*cI.^2 - 1); % [rad/s]

    % Unit conversions
    rad2deg = 180/pi;
    sec2day = 86400;

    rates.Omega_rad_s = dOmega;
    rates.omega_rad_s = domega;
    rates.M_rad_s     = dM;

    rates.Omega_deg_day = dOmega * rad2deg * sec2day;
    rates.omega_deg_day = domega * rad2deg * sec2day;
    rates.M_deg_day     = dM     * rad2deg * sec2day;
end
