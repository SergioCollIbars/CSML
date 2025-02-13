function [planetParams, poleParams, Cmat_true, Smat_true, ...
    TIME, DOM] = load_universe(system, tspan, frec)
    %%                    LOAD UNIVERSE FUNCTION
    % Description: Based on the especified system, load the planet, pole 
    % params. Also load the initial conditions.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024

    tmin = tspan(1);
    tmax = tspan(2);

    planetParams = zeros(1, 9);
    poleParams   = zeros(1, 8);

    R_earth = cspice_bodvrd('EARTH', 'RADII', 3); % [Km]
    R_moon  = cspice_bodvrd('MOON', 'RADII', 3);  % [Km]

    GM_earth = cspice_bodvrd('EARTH', 'GM', 1);   % [km^3/s^2]
    GM_moon  = cspice_bodvrd('MOON', 'GM', 1);    % [km^3/s^2]

    GM1 = GM_earth * 1E9;    % [m^3/s^2]
    GM2 = GM_moon  * 1E9;    % [m^3/s^2]
    
    path1 = "HARMCOEFS_EARTH_1.txt";
    path2 = "HARMCOEFS_MOON_1.txt";

    [Cmat1, Smat1, ~] = readCoeff(path1); % grav. field primary
    [Cmat2, Smat2, ~] = readCoeff(path2); % grav. field secondary

    Cmat_true = {Cmat1, Cmat2};
    Smat_true = {Smat1, Smat2};

    normalized = 1;
    if(system == "CR3BP"), n_max = 0; else, n_max = 8; end

    % system paramters. Earth & Moon
    D = 384399e3;  % [m]
    planetParams(1) = GM2 / (GM1 + GM2); % mass ratio
    planetParams(2) = D;                 % primaries distance

    % define time dimensionalization
    n = sqrt((GM1 + GM2) / D^3); % mean motion circular orbit [1/s]
    planetParams(3) = n;
    planetParams(4) = R_earth(1) * 1E3;       % reference radius primary   
    planetParams(5) = R_moon(1)  * 1E3;       % reference radius secondary
    planetParams(6) = n_max;     % max SH zonal
    planetParams(7) = normalized;% normalized grav. coeff. option  
    planetParams(8) = GM_earth * 1E9;       % primary point mass parameter [m^3/s^2]
    planetParams(9) = GM_moon  * 1E9;       % secondary point mass parameter [m^3/s^2]
    planetParams(10) = 1.3;                 % SRP scaling factor
    planetParams(11) = 1000;                % S/C mass [Kg]
    planetParams(12) = 50 / (D^2);          % S/C area [-]

    % define pole parameters w.r.t J2000. Earth and Moon
    poleParams(1) =  -pi/2;                   % RA Earth
    poleParams(2) = deg2rad(90);              % DEC Earth
    poleParams(3) = 0;                        % prime meridian Earth
    poleParams(4) = 0;                        % ang. velocity Earth
    poleParams(5) =  -pi/2;                   % RA Earth
    poleParams(6) = deg2rad(90);              % DEC Moon
    poleParams(7) = 0;                        % prime meridian Moon
    poleParams(8) = 0;                        % ang. velocity Moon
    
    % time vector [-]
    TIME = [tmin, tmax];
    if(system == "EPHEM")
        ET = cspice_str2et({'2020-01-10 12:00:00 TDB'});
        TIME  = TIME  + ET(1)*planetParams(3);  % [-]
    end
    
     % Date Of Measurement vector
     n = round((TIME(end) - TIME(1))*(frec/planetParams(3)) + 1);
     DOM = linspace(TIME(1), TIME(end), n);

end

