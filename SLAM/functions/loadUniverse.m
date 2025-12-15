function [planetParams,poleParams, Cnm, Snm, ...
    initCond, t_range] = loadUniverse(metaData_file)
    % Description: Load the planet / asteorid information and the
    % orientation information for the body
    mtd = readParams("data/"+metaData_file);
    params = readParams("data/"+mtd.folder+"/Universe.txt");
    
    % Bennu parameters
    GM = params.GM;
    n_max = params.n_max;

    path = params.path;
    [Cnm, Snm, Re] = readCoeff(path);
    normalized = params.normalized;
    
    % orientation parameters
    W = params.W;  % Rotation ang. vel   [rad/s]
    W0 = params.W0;                   % Initial asteroid longitude
    RA = deg2rad(params.RA);    % Right Ascension     [rad]
    DEC = deg2rad(params.DEC);  % Declination         [rad]

    poleParams   = [W, W0, RA, DEC];
    planetParams = [GM, Re, n_max, normalized];

    % initial conditions
    r      = params.r;      % [m]
    phi    = deg2rad(params.phi);
    lambda = params.lambda;
    theta  = pi/2 - phi;% Orbit colatitude [m]
    R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
        sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
        cos(theta), -sin(theta), 0];
    r0 = R * [r;0;0];           % [ACI]
    v0 = R * [0;0;sqrt(GM/r)];  % [ACI]
    
    initCond = [r0;v0];

    % Initial and final time
    t_range = [params.tmin*86400, params.tmax*86400];
end

