function [AccObj] = acc_module(AttitudeObj, OrbitObj, AccObj, Time,  save,...
                    name, readData)
    %%                      ACCELERATION MODULE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This module is in charge of compute the accelerometer
    %   measurements along the orbit. 
    %
    %   Input:
    %       OrbitObj:  orbit class object
    %       AttitudeObj: attitude class object
    %       AccObj: acceleration class object
    %       Time: time vector
    %       save: save variables option boolean
    %       name: saved file name 
    %
    %   Output:
    %       AccObj: acceleration class object
    % --------------------------------------------------------------------%
    disp('  Acceleration module')

    if(readData == 1)
        % parse data
        AccObj = parseDataAcc(AccObj);

        ri = AccObj.ri;
        omega = AccObj.omega;
        Omega = AccObj.Omega;
        Rtot = AccObj.Rtot;
        
        % select C / S and W varaibles
        OrbitObj = OrbitObj.C_Wmatrix();
    else
        % obtain orbit variables
        ri = OrbitObj.ri;                                      % [m]
    
        % obtain attitude variables
        omega = AttitudeObj.omega;                             % [rad / s]
        Omega = AttitudeObj.Omega;                             % [rad / s^2]
       
        % rotation values
        Rtot = OrbitObj.BN;
    end
    
   
    % Orbit data
    GM = OrbitObj.OrbitData(1);                              % [m^3 / s^2]
    Re = OrbitObj.OrbitData(10);
    n_max = OrbitObj.OrbitData(9);
    Normalized = OrbitObj.OrbitData(13);

    % Orbit and Planet parameters
    C_mat = OrbitObj.C;
    S_mat = OrbitObj.S;
    
    % Rotation parameters
    W = OrbitObj.W;
    W0 = OrbitObj.W0;
    RA = OrbitObj.RA;
    DEC = OrbitObj.DEC;

    % time steps number
    T = length(ri);

    % compute gravity gradient in body axis    
    dUb = zeros(3, T);
    
    % compute gravity tensor in body axis
    accT = zeros(3*T, 3);

    % bias value. Truth
    B = [-282,8.88,16.42;1848,-636,-19500;5,9699,-113] ...
        * 1E-9/2;  % [1/s^2]
    

    % drift value. Truth
    D = [0.0065,0.0103, 0.0001; 3.123,0.195,0.9374;0.001, 2.98,0.0033] ...
        * 1E-9/(2*86400); % [1/s^2/s]
    
%     B = zeros(3, 3);
%     D = B;

    % noise value
    noise = zeros(9, T);
    At = Time(2) - Time(1);
    Fs = 1/(At);
    sigma_ii = AccObj.sigmaN_ii;
    sigma_ij = AccObj.sigmaN_ij;
    if(AccObj.noiseProfile == 1)
        noise = randn(9, T);
        noise = noise.*[sigma_ii;sigma_ij;sigma_ij;sigma_ij;...
            sigma_ii;sigma_ij;sigma_ij;sigma_ij;sigma_ii];
    elseif(AccObj.noiseProfile == 2)
        for i = 1:9
            if(i == 1 || i == 5 || i == 9)
                wn = sigma_ii;
                s = 6.32E-12 / wn;      % scaling value
                pn = wn / (50 * s);
            else
                wn = sigma_ij;
                s = 6.32E-12 / wn;      % scaling value
                pn = wn / (50 * s);
            end
            [noise(i, :), ~, ~] = noise_profile(wn, pn, T, Fs);
        end
    else
        warning('noise profile value not understand')
    end

    tic
    for j =1:T
        % Body rotation matrix @ current time step
        Rj = Rtot(3*j-2 : 3*j, :);
        
        % ECEF rotation @ current time step
        Wt = W * Time(j) + W0;
        BN =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        NB = BN';

        % compute ang vel and acc matrices. Body frame
        [omega2M, OmegaM] = angularMatrix(omega(:, j), Omega(:, j));

        % CoM postion ECEF
        r_ECEF = BN * ri(:, j);
        
        % compute gravity quantities. ECEF
        [~, dU_ECEF, ddU_ECEF] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                r_ECEF, Re, GM, ...
                                                Normalized);
        % rotate 2 ECI
        ddU_ECI = NB * ddU_ECEF * BN;
        dU_ECI = NB * dU_ECEF;

        % rotate 2 body 
        ddU = Rj * ddU_ECI * Rj';
        dU = Rj * dU_ECI;

        % include 3BP gravity tensor. Hill approximation
        if(OrbitObj.OrbitData(12) == 1)
            F = sqrt(OrbitObj.GMs / OrbitObj.alphaK(1)^3);
            ddU_3BP = [3*F*F, 0, 0;0 , 0, 0;0, 0, -F*F];
            SAR_N = OrbitObj.SAR_N(3*j-2:3*j, :);
            ddU_3BP = SAR_N' * ddU_3BP * SAR_N;
        else
            ddU_3BP = zeros(3,3);
        end

        % compute acc tensor
        up = 3*j;
        down = up - 2;
        accT(down:up, :) = -ddU + omega2M + OmegaM - ddU_3BP;

        % save potential body frame
        dUb(:, j) = dU;
        
        % bias term
        if( j ~= 1)
            B = B + D * At;
        end

        % generate random noise
        accT(down:up, :) = accT(down:up, :) + ...
                             [noise(1, j), noise(2, j), noise(3, j);...
                             noise(4, j), noise(5, j), noise(6, j);...
                             noise(7, j), noise(8, j), noise(9, j)] + ...
                             B;
    end
    toc

    % save values into object
    AccObj.accT = accT;
    AccObj.dUb = dUb;
    AccObj.t = Time;

    % save data in file 
    if(save == true)
        saveData_acc(AccObj, name);
    end
end