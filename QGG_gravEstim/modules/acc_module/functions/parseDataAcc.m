function [AccObj] = parseDataAcc(AccObj)
    %%                       PARSE DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 23/02/2023
    %
    %   Description: This function loads and parse data from the data files
    %   generated. 
    %
    %   Input:
    %       AccObj:  acceleration class object
    %
    %   Output:
    %       AccObj:  acceleration class object
    % --------------------------------------------------------------------%
    
    % Inputs
    addpath('./data_files/');

    % read data files    
    name = "orbitData.txt";
    T_orbit = readtable(name);
    
    name = "attitudeData.txt";
    T_att = readtable(name);

    name = "ECI2Body.txt";
    T_BN = readtable(name);

    t = T_orbit.TIME;

    % create data matrices
    ri = zeros(3, length(t));                    % Inertial pos CoM
    vi = zeros(3, length(t));                    % Inertial vel CoM
    
    angVel = zeros(3, length(t));                % Ang vel Body
    angAcc = zeros(3, length(t));                % Ang acc Body

    BN = zeros(3*length(t), 3);                  % Rotation matrix

    % parse data
    ri(1, :) = T_orbit.ri_x;                            
    ri(2, :) = T_orbit.ri_y;
    ri(3, :) = T_orbit.ri_z;
    
    vi(1, :) = T_orbit.vi_x;
    vi(2, :) = T_orbit.vi_y;
    vi(3, :) = T_orbit.vi_z;
    
    angVel(1, :) = T_att.omega_x;
    angVel(2, :) = T_att.omega_y;
    angVel(3, :) = T_att.omega_z;
    
    angAcc(1, :) = T_att.Omega_x;
    angAcc(2, :) = T_att.Omega_y;
    angAcc(3, :) = T_att.Omega_z;

    for k = 1:length(t)
        % compute tensor matrix
        up = 3 * k;
        down = up - 2;

        % compute rotation matrix
        BN(down:up, :) = [T_BN.BN_1(down:up), T_BN.BN_2(down:up), ...
            T_BN.BN_3(down:up) ];
    end

    % load to estimation class
    AccObj.ri = ri;
    AccObj.vi = vi;

    AccObj.omega = angVel;
    AccObj.Omega = angAcc;

    AccObj.Rtot = BN;
end

