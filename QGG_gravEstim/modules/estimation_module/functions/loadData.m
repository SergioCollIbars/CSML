function [GG_meas, pos_meas, vel_meas, angVel_meas, ...
    angAcc_meas, B_ACI, TIME] = loadData()
    %%                       LOAD DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 23/02/2023
    %
    %   Description: This function loads and parse data from the data files
    %   generated. 
    %
    %   Input:
    %       estimObj:  estimation class object
    %
    %   Output:
    %       estimObj:  estimation class object
    % --------------------------------------------------------------------%
    
    % Inputs
    addpath('./data_files/');

    % read data files
    name = "accData.txt";
    T_acc = readtable(name);
    
    name = "orbitData.txt";
    T_orbit = readtable(name);
    
    name = "attitudeData.txt";
    T_att = readtable(name);

    name = "ECI2Body.txt";
    T_BN = readtable(name);

    t = T_acc.TIME;

    % create data matrices
    accT = zeros(3*length(t), 3);        % Acc tensor Body frame
    
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
        accT(down, 1) = T_acc.ad_xx(k);
        accT(down, 2) = T_acc.ad_xy(k);
        accT(down, 3) = T_acc.ad_xz(k);

        accT(down + 1, 1) = T_acc.ad_yx(k);
        accT(down + 1, 2) = T_acc.ad_yy(k);
        accT(down + 1, 3) = T_acc.ad_yz(k);

        accT(down + 2, 1) = T_acc.ad_zx(k);
        accT(down + 2, 2) = T_acc.ad_zy(k);
        accT(down + 2, 3) = T_acc.ad_zz(k);

        % compute rotation matrix. Body to ACI
        BN(down:up, :) = [T_BN.BN_1(down:up), T_BN.BN_2(down:up), ...
            T_BN.BN_3(down:up) ];
    end
    
    % save measurements
    GG_meas = accT;
    pos_meas = ri;
    vel_meas = vi;
    angVel_meas = angVel;
    angAcc_meas = angAcc;
    B_ACI = BN;
    TIME = t;
end

