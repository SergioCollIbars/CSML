function [planetParams,poleParams, Cnm, Snm, ...
    initCond] = loadUniverse()
    % Description: Load the planet / asteorid information and the
    % orientation information for the body
    
    % Bennu parameters
    GM = 5.2;
    n_max = 6;

    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    normalized = 1;
    
    % orientation parameters
    W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
    W0 = 0;                   % Initial asteroid longitude
    RA = deg2rad(86.6388);    % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);  % Declination         [rad]

    poleParams   = [W, W0, RA, DEC];
    planetParams = [GM, Re, n_max, normalized];

    % initial conditions
    r      = 0.7E3;      % [m]
    phi    = pi/2;
    lambda = 0;
    theta  = pi/2 - phi;% Orbit colatitude [m]
    R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
        sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
        cos(theta), -sin(theta), 0];
    r0 = R * [r;0;0];           % [ACI]
    v0 = R * [0;0;sqrt(GM/r)];  % [ACI]
    
    initCond = [r0;v0];
end

