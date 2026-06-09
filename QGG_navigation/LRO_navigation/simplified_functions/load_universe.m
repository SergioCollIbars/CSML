function [planetParams, Cmat_true, Smat_true] = load_universe()
    %%                    LOAD UNIVERSE FUNCTION
    % Description: Based on the especified system, load the planet, pole 
    % params. Also load the initial conditions.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024

    R      = cspice_bodvrd('EARTH', 'RADII', 3); % [Km]
    R_E    = R(1)*1e3; % [m]
    R_M    = 1738*1E3; % [m]

    GM_earth = cspice_bodvrd('EARTH', 'GM', 1);   % [km^3/s^2]

    GM_E   = GM_earth * 1E9;    % [m^3/s^2]
    GM_M   = 4.9028001224453001E12;

    [GM_S] = cspice_bodvrd('SUN', 'GM', 1)*1E9; 
    [R]    = cspice_bodvrd('SUN', 'RADII', 3).*1E3;
    R_S    = R(1);
    
    path1 = "HARMCOEFS_EARTH_1.txt";
    path2 = "HARMCOEFS_MOON_GRGM1200.txt";

    [Cmat1, Smat1, ~] = readCoeff(path1); % grav. field primary
    [Cmat2, Smat2, ~] = readCoeff(path2); % grav. field secondary

    Cmat_true = {Cmat1, Cmat2};
    Smat_true = {Smat1, Smat2};

    normalized = 1;
    n_maxE = 0; 
    n_maxM = 400; 

    planetParams = [GM_E, GM_M, R_E, R_M, n_maxM, normalized,...
        GM_S,n_maxE,R_S];
end

