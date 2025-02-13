function [OrbitObj, AttitudeObj] = orbit_module(OrbitObj, AttitudeObj, t,...
    save)
    %%                          ORBIT MODULE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This module is in charge of compute the orbital
    %   movement for a prescribed orbit. 
    %
    %   Input:
    %       OrbitObj:  orbit class object
    %       t: time vector
    %       save: save variables option boolean
    %
    %   Output:
    %       OrbitObj:  orbit class object
    % --------------------------------------------------------------------%
    disp('  Orbit module')
    % Ideal rotation?
    ideal = false;
    if (ideal), disp('      SC ideal orbit rotation on'); end

    % gravity parameter
    mu = OrbitObj.OrbitData(1);                      % [m^3 / s^2]

    % orbital parameters
    e = OrbitObj.OrbitData(2);                       % eccentrycity
    a = OrbitObj.OrbitData(3);                       % semi major axis [m]
    rho = a * (1 - e^2);                             % orbital param [m]

    i = deg2rad(OrbitObj.OrbitData(4));              % inclination [rad]
    omega = deg2rad(OrbitObj.OrbitData(5));          % arg periapsis [rad]
    Omega = deg2rad(OrbitObj.OrbitData(6));          % RAAN [rad]
    f = deg2rad(OrbitObj.OrbitData(7));              % true anomaly [rad]
    
    % compute initial state vectors (orbit frame)
    r0 = rho / (1 + e * cos(f)) * [cos(f);...
                                  sin(f);...
                                  0];
    v0 = sqrt(mu / rho) * [-sin(f);...
                          e + cos(f);...
                          0];
    
    % compute rotation matrix: Body to ACI
    [BN] = rotationMatrix(Omega, i, omega, [3,1,3]);

    % rotate initial state vectors. {I, J, K} frame
    r0 = BN' * r0;
    v0 = BN' * v0;
    
    % define initial conditions
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];

    % define integration options
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    
    % Select harmonic coefficients matrix and rotation params
    OrbitObj = OrbitObj.C_Wmatrix();
    
    % planet pole parameters
    W = OrbitObj.W;
    W0 = OrbitObj.W0;
    RA = OrbitObj.RA;
    DEC = OrbitObj.DEC;

    % ODE 113
    [~, state] = ode113(@(t, x) EoM(t, x, mu, OrbitObj), t, X0, options);
    state = state';
    
    % body frame matrix definition
    rb = zeros(3, length(t));
    vb = zeros(3, length(t));

    % orbit elements matrix definition
    alpha = zeros(8, length(t));
    ang = zeros(3, length(t));
    ang_dot = zeros(3, length(t));
    
    % total rotation matrix
    BN = zeros(3 * length(t), 3);
    ACAF_N = zeros(3 * length(t), 3);
    SAR_N = zeros(3 * length(t), 3);

    % Sun position definition
    rs_ACI = zeros(3, length(t));

    tic
    for k = 1:length(t)
        % inertial state vectors
        ri = state(1:3, k);
        vi = state(4:6, k);
        
        % compute orbital elements
        alpha(:, k) = orbitalElem(ri, vi, mu);

        % change orbit sense rotation
        if(isnan(alpha(7, k)) || isnan(alpha(8, k)) || OrbitObj.OrbitData(9) < 2)
           BN1 = rotationMatrix(Omega, i, omega, [3,1,3]);
        else
           BN1 = rotationMatrix(alpha(7, k), alpha(6, k), alpha(8, k), [3,1,3]);
        end

        % compute orbit frame
        ro = BN1 * ri;
        vo = BN1 * vi;
        
        % compute ideal or real SC rotation matrix
        if(ideal == 0)
            BN2 = rotationMatrix(AttitudeObj.theta1(k) + f, ...
                AttitudeObj.theta2(k), AttitudeObj.theta3(k), [3,1,3]);
            identM = eye(3, 3);

            if(isequal(BN2, identM))
                BN1 = eye(3, 3);
            end
        else
            theta_dot = cross(ro, vo) / (vecnorm(ro)^2);
            theta_dot_u = theta_dot / vecnorm(theta_dot); 

            h = cross(ro, vo);
            theta_ddot2  = -2 * vecnorm(h) / (vecnorm(ro)^4) * ...
                dot(ro, vo) * theta_dot_u;
 
            theta = atan2(ro(2), ro(1));

            BN2 = rotationMatrix(theta, 0, 0, [3,1,3]);

            ang(:, k) = theta_dot;
            ang_dot(:, k) = theta_ddot2;
        end
        
        % compute body frame state vectors
        rb(:, k) = BN2 * BN1 * ri;
        vb(:, k) = BN2 * BN1 * vi;

        % compute rotation matrix
        BN(3*k-2:3*k, :) = BN2 * BN1;
        
        % Inertial to ACAF rotation 
        Wt = W * t(k) + W0;
        R = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        ACAF_N(3*k-2:3*k, :) = R;

        % Inertial to SAR rotation
        [r_ACI, v_ACI] = Kepler_motion(OrbitObj.alphaK, ...
            OrbitObj.GMs, t(k));
        H = cross(r_ACI, v_ACI);
        H = H / vecnorm(H);
    
        D = r_ACI./vecnorm(r_ACI);
        F = cross(H, D);

        SAR_N(3*k-2:3*k, :) = [D';F';H'];
        rs_ACI(:, k) = r_ACI;

    end
    toc

    if(ideal == 1)
        AttitudeObj.omega = ang;
        AttitudeObj.Omega = ang_dot;
    end

    % save variable into orbit object
    OrbitObj.ri = state(1:3, :);
    OrbitObj.vi = state(4:6, :);

    OrbitObj.rb = rb;
    OrbitObj.vb = vb;
    OrbitObj.rs_ACI = rs_ACI;
    OrbitObj.t = t;
    
    OrbitObj.alpha = alpha;
    OrbitObj.BN = BN;
    OrbitObj.ACAF_N = ACAF_N;
    OrbitObj.SAR_N = SAR_N;
    
    % save data in file 
    if(save == true)
        saveData_orbit(OrbitObj, "orbitData.txt");
        saveData_attitude(AttitudeObj, "attitudeData.txt")
    end

end