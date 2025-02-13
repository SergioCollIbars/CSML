function [X0] = load_initCond(system, planetParams, TIME)
    %%                    LOAD INITIAL CONDITIONS FUNCTION
    % Description: Based on the especified system, load the initial
    % conditions of the orbit.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024

    if(system == "2BP") % point mass 2BP dynamics. Earth orbit
        R = planetParams(2);
        GM = planetParams(1);

        % Initital conditions. Orbital parameters
        e = 0.6;                     % eccentrycity
        a = R + 10000e3;               % semi major axis [m]
        rho = a * (1 - e^2);         % orbital param [m]
        
        i = deg2rad(45);              % inclination [rad]
        omega = deg2rad(0);          % arg periapsis [rad]
        Omega = deg2rad(0);          % RAAN [rad]
        f = deg2rad(-180);           % true anomaly [rad]
        
        [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
    elseif(system == "F2BP")
        GM = planetParams(1);

        e = 0;                       % eccentrycity
        a = 2000;                    % semi major axis [m]
        rho = a * (1 - e^2);         % orbital param [m]
        
        i = deg2rad(45);             % inclination [rad]
        omega = deg2rad(0);          % arg periapsis [rad]
        Omega = deg2rad(0);          % RAAN [rad]
        f = deg2rad(-180);              % true anomaly [rad]

        [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
    elseif(system == "CR3BP" || system == "FCR3BP")
         r0 = [1.021968177072928; 0; -0.18206];
         v0 = [0; -0.1031401430288178; 0]; % L1 orbit

         [r0, v0] = rotate2inertial(r0, v0, 0, 1);  % [-] and [-]
    elseif(system == "EPHEM")
        et = TIME(1) / planetParams(3);
        [state, ~] = cspice_spkezr('-60000', et, 'J2000', 'NONE', '3');

        r0 = state(1:3);
        v0 = state(4:6);

        r0 = r0*1E3/planetParams(2);
        v0 = v0*1E3/(planetParams(3)*planetParams(2));
    end

    % define initial conditions
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
end

