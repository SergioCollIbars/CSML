function [planetParams, Cnm_list, Snm_list, ...
    initCond, t_range] = loadUniverse(folder_Name)
    % Description: Load the planet / asteorid information and the
    % orientation information for the body
    params = readParams("data/"+folder_Name+"/Universe.txt");
    
    % Bennu parameters
    n_max = params.n_max;
    normalized = params.normalized;

    % Get GM for the Moon [m^3/s^2]
    [GM_M] = cspice_bodvrd('MOON', 'GM', 1)*1E9; 

    % get Moon Radius [m]
    R_M  = 1738*1E3;
    
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

    % Initial and final time [sec]
    et        = cspice_str2et(params.utc_start);
    t_init    = et;
    t_final   = et + params.tmax*86400; 
    t_range   = [t_init, t_final];

    % initial conditions
    tgt       = params.tgt;
    if(tgt == "LRO")
        tgt       = 'LUNAR RECONNAISSANCE ORBITER';
        observer  = 'MOON';
        ref_frame = 'J2000';
        [sc_SPICE, ~] = cspice_spkezr(tgt, et, ...
                        ref_frame, 'NONE', observer);

        initCond = sc_SPICE.*1E3;   % [m] and [m/s]
    elseif(tgt == "CUSTOM")
        % TBD. Custom trajectory
        a    = params.a*1E3;            % [m]
        e    = params.e;                % [-]
        i    = deg2rad(params.i);       % [rad]
        RAAN = deg2rad(params.RAAN);    % [rad]
        omega= deg2rad(params.omega);   % [rad]
        f    =deg2rad(params.f);        % [rad] 
        
        rho  = a * (1 - e^2);           % [m] 
        
        % Initial conditions Orbit Frame (OF)
        r0_OF = rho / (1 + e * cos(f)) * [cos(f);sin(f);0];
        v0_OF = sqrt(GM_M / rho) * [-sin(f);e + cos(f);0];
        
        % compute rotation matrix: OF to Inertial
        [BN] = rotationMatrix(RAAN, i, omega, [3,1,3]);
    
        % rotate initial state vectors.
        r0 = BN' * r0_OF;
        v0 = BN' * v0_OF;

        initCond = [r0;v0];     % [m] and [m/s]
    else
        error('Orbit type not defined');
    end
end

