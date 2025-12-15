function [planetParams, Cnm_list, Snm_list, ...
    initCond, t_range] = loadUniverse(metaData_file)
    % Description: Load the planet / asteorid information and the
    % orientation information for the body
    mtd = readParams("data/"+metaData_file);
    params = readParams("data/"+mtd.folder+"/Universe.txt");
    
    % Bennu parameters
    n_max = params.n_max;
    normalized = params.normalized;

    % Get GM for the Moon [m^3/s^2]
    [GM_M] = cspice_bodvrd('MOON', 'GM', 1)*1E9; 

    % get Moon Radius [m]
    [R]  = cspice_bodvrd('MOON', 'RADII', 3)*1E3;
    R_M  = R(1);
    
    % Get GM for the Earth [m^3/s^2]
    [GM_E] = cspice_bodvrd('EARTH', 'GM', 1)*1E9; 

    % get Earth Radius [m]
    [R]  = cspice_bodvrd('EARTH', 'RADII', 3)*1E3;
    R_E  = R(1);

    % get gravity coefficients
    path_Earth = params.path_Earth;
    [Cnm_E, Snm_E, ~] = readCoeff(path_Earth);

    path_Moon = params.path_Moon;
    [Cnm_M, Snm_M, ~] = readCoeff(path_Moon);

    % save in list
    Cnm_list = {Cnm_E, Cnm_M};
    Snm_list = {Snm_E, Snm_M};

    planetParams = [GM_E, GM_M, R_E, R_M, n_max, normalized];

    % initial conditions
    tgt       = params.tgt;
    observer  = params.observer;
    ref_frame = params.ref_frame;
    et        = cspice_str2et(params.utc_start);
    [sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);
    
    initCond = sc_SPICE.*1E3;   % [m] and [m/s]

    % Initial and final time [sec]
    t_init  = et + params.tmin*86400;
    t_final = et + params.tmax*86400; 
    t_range = [t_init, t_final];
end

