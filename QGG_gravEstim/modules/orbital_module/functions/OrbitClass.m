classdef OrbitClass
    %%                      ORBIT CLASS
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: Class to store all the information related with the
    %   orbit module
    %
    % --------------------------------------------------------------------%
    properties
        % profile configuration
        OrbitData
        SH_path
    
        % state variables
        ri                          % inertial frame position
        vi                          % inertial frame velocity
        rb                          % body frame position
        vb                          % body frame velocity
        rs_ACI                      % Sun position. ACI frame
        alpha                       % orbital elements vector
        alphaK                      % orbital elements keplerian motion

        % time vector
        t
        
        % rotation parameters
        W_Earth = 0.7292115E-4;    % Earth angular velocity [rad / s]
        W0_Earth = 0;              % Earth prime meridian at t0 [rad]

        W_Bennu = 4.06130329511851E-4; % Bennu angular velocity [rad / s]
        W0_Bennu = 0;              % Bennu prime meridian at t0 [rad]

        W_Eros = 0.000331165979836797; % Eros angular velocity [rad / s]
        W0_Eros = 0;              % Eros prime meridian at t0 [rad]

        % rotation matrices
        BN                         % Inertial to Body frame
        SAR_N                      % Inertial to SAR frame
        ACAF_N                     % Inertial to ACAF frame

        % Cnm matrix. Perturbed coefficients
        Cnm_Earth = [1, 0, 0, 0;...
                     0, 0, 0, 0;...
                    -1.08E-3, 0, 1.57E-6, 0;...
                     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 

        Snm_Earth = [0, 0, 0, 0;...
                     0, 0, 0, 0;...
                     0, 0, -9.03E-7, 0;...
                     0, 2.68E-7, -2.12E-7, 1.98E-7]; 


        % C & S. Harmonic coefficients selected
        C
        S

        % W & W0 selected values
        W
        W0

        % Sun model
        Rs = 6.96E8                 % Sun radius [m]
        Ts = 5778;                  % Sun temperature [k]
        rcE = 1.495978707E11;       % distance Earth to Sun
        rcB = 1.1259 * ...
            1.495978707E11;         % distance Bennu to Sun
        rcEr = 1.45 * ...
            1.495978707E11;         % distance Eros to Sun
        rc                          % distance form planet to Sun
        GMs = 1.32712440018E20;     % Gravity param Sun

        % Planet Keplerian motion (a, e, i, RAAN, omega, f)
        alphaBennu = [1.1259*1.495978707E11, 0.20372, ...
            deg2rad(6.034),deg2rad(2.017), deg2rad(66.304), ...
            deg2rad(64.541)];

        alphaEros = [1.45*1.495978707E11, 0.22, ...
            deg2rad(10.8),deg2rad(0), deg2rad(0), ...
            deg2rad(0)];

        alphaEarth = [1.00000011*1.495978707E11, 0.0167102, ...
            deg2rad(0.00005),deg2rad(-11.26064), deg2rad(102.94719), ...
            deg2rad( 100.46435)];

        % ACAF parameters
        RA;                      % Right ascension [rad]
        DEC;                     % declination [rad]


        % physics constants
        c = 2.99792458E8;           % speed of light [m / s]
        gamma = 5.67E-8;            % Stefan-Boltzman constant [w*m^-2*k^4]
    end
    methods
        function obj = C_Wmatrix(obj)
                % Cnm length
                if(obj.OrbitData(8) == 2)
                    obj = readCoeff(obj);

                    obj.W = obj.W_Bennu;
                    obj.W0 = obj.W0_Bennu;

                    obj.rc = obj.rcB;

                    obj.alphaK = obj.alphaBennu;

                    obj.RA  = deg2rad(86.6388);
                    obj.DEC = deg2rad(-65.1086);

                elseif(obj.OrbitData(8) == 1)
                    obj = readCoeff(obj);

                    obj.W = obj.W_Earth;
                    obj.W0 = obj.W0_Earth;

                    obj.rc = obj.rcE;

                    obj.alphaK = obj.alphaEarth;
               elseif(obj.OrbitData(8) == 3)
                    obj = readCoeff(obj);

                    obj.W = obj.W_Eros;
                    obj.W0 = obj.W0_Eros;

                    obj.rc = obj.rcEr;

                    obj.alphaK = obj.alphaEros;

                    obj.RA  = deg2rad(11.363);
                    obj.DEC = deg2rad(17.232);
                end
        end

        function obj = readCoeff(obj)
            list = table2array(readtable(obj.SH_path));
            degree = list(1);
            X = list(4:end);
            obj.OrbitData(13) = list(3);
            
            % count number of coefficients
            Nc = 1;
            for k = 2:degree
                Nc = Nc + k + 1;
            end
            Ns = 0;
            for k = 2:degree
                Ns = Ns + k;
            end
            
            % define matrices
            Cmat = zeros(degree + 1, degree + 1);
            Smat = zeros(degree + 1, degree + 1);
            Cmat(1, 1) = 1;

            n = 2;
            m = 0;
            for j = 2:Nc
                N = n + 1;
                M = m + 1;
       
                Cmat(N, M) = X(j);
                if(m < n)
                    m = m + 1;
                else
                    m = 0;
                    n = n +1;
                end
            end
        
            n = 2;
            m = 0;
            for j = Nc + 1:Ns + Nc
                N = n + 1;
                M = m + 2;
        
                Smat(N, M) = X(j);
                if(m < n - 1)
                    m = m + 1;
                else
                    m = 0;
                    n = n + 1;
                end
           
            end
            % assign to C and S matrices
            obj.C = Cmat;
            obj.S = Smat;
        end
    end
end