classdef estimationClass
    %%                      ESTIMATION CLASS
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 23/02/2023
    %
    %   Description: Class to store all the information related with the
    %   estimation module
    %
    % --------------------------------------------------------------------%
    %% PROPERTIES
    properties
        % profile configuration
        EstimData

        % CONFIG VALUES
        GM                      % Gravity constant
        Re                      % planet radious
        poleParams              % planet pole parameters
        C_mat                   % C harmonics coefficients
        S_mat                   % S harmonics coefficients
        X0                      % init nominal state
        delta_X0                % init deviation
        atm                     % atmosphere params
        RS                      % position to Sun
        
        % MEASUREMENTS
        meas                    % measurements vector until DCO
        measT                   % total meas. 
        TIME                    % time vector
        t0                      % init time at the epoch
        extraParam              % extra parameter mode

        % INIT VARIABLES
        n_max                   % Maximum SH degree
        Nc                      % Number of Cnm states.
        Ns                      % Number of Snm states.
        Np                      % Number of params 
        Nt                      % Number of time steps
        Nm                      % Number of meas per time step

        Xfilter                 % filter state estimation
        Xfilterj                % filter state estimation per iter
        pert                    % estimate perturbation
        deltaXj                 % estimate state peturb per iter
        deltaX
        sigma_t                 % STD over time 

        Hi                      % visibility mat @ time step
        Pj                      % filter covariance per iter
        Pj0                     % a posteriori covariance at epoch
        QTilde0                 % inital process noise continous
        Bw                      % perturbed acc time matrix

        prefit                  % prefit vector
        postfit                 % postfit vector
        prefitj                 % prefit vector per iter
        postfitj                 % posfit vector per iter

        P0                      % initial state covariance
        P0_inv                  % initial state covariance inverse

        R0                      % inital meas covariance
        R0_inv                  % initial meas covariance inverse

        Xnot                    % acumulated deviation
    end
    %% METHODS
    methods
        function obj = initVariables(obj, TIME, MaxIter)
            % fill variables
            obj.TIME = TIME;
            obj.Np = length(obj.X0);
            obj.Nt = length(obj.TIME);
            obj.Nm = 6;

            obj.Xfilter = zeros(obj.Np, obj.Nt);
            obj.Xfilterj = ones(obj.Np*MaxIter, obj.Nt) * NaN;
            obj.pert = zeros(obj.Np, obj.Nt);
            obj.deltaXj = ones(obj.Np*MaxIter, obj.Nt) * NaN;

            obj.Hi =  zeros(obj.Nm*obj.Nt, obj.Np);
            obj.Pj = zeros(obj.Np*obj.Nt, obj.Np, MaxIter);
            obj.Pj0 = zeros(obj.Np, obj.Np, MaxIter);
            obj.sigma_t = zeros(obj.Np, obj.Nt);

            obj.prefit = ones(obj.Nm, obj.Nt)*NaN;
            obj.postfit = ones(obj.Nm, obj.Nt)*NaN;

            obj.prefitj = ones(obj.Nm*MaxIter, obj.Nt) * NaN;
            obj.postfitj = ones(obj.Nm*MaxIter, obj.Nt) * NaN;

            obj.Xnot = obj.delta_X0;
            obj.deltaX = obj.delta_X0;
            
            % initialize covariances
            sigma2_Cnm = ones(1, obj.Nc)*(1)^2;
            sigma2_Cnm(1) = obj.EstimData(7)^2;
            sigma2_Snm = ones(1, obj.Ns)*(1)^2;
            obj.P0 = diag([sigma2_Cnm, sigma2_Snm]);

            if(obj.extraParam ~= 0)
% %                 % phase A
% %                 sigma2_E = [(1E-10)^2, (1E-15)^2];
% %                 sigma2_E = repmat(sigma2_E, [1, 6]);
% %                 b = 1E-9;
% %                 sigma2_E = [b, 1E-16, b, 1e-16, b, 1E-19, b, 1E-15,...
% %                     b, 1E-15, b, 1E-16].^2;

                % phase B
                sigma2_E = [3.63432250269746e-14,...
                            5.20404953330448e-19,...
                            1.44897421768161e-12,...
                            2.07486466975169e-17,...
                            1.44897421767654e-11,...
                            2.07486466975336e-19,...
                            3.63432250835262e-9,...
                            5.20404953474790e-20,...
                            1.44897421767758e-20,...
                            2.07486466975374e-16,...
                            3.63432251628775e-20,...
                            5.20404953318338e-17].^2;

                
                obj.P0 = diag([sigma2_Cnm, sigma2_Snm, sigma2_E]);
            end

            obj.P0_inv = diag([diag(1./obj.P0)]);
            
            sigma2_ii = obj.EstimData(3)^2;
            sigma2_ij = obj.EstimData(4)^2;
            obj.R0 = diag([sigma2_ii, sigma2_ij, sigma2_ij,...
                sigma2_ii, sigma2_ij, sigma2_ii]);
            obj.R0_inv = inv(obj.R0);

            % process noise covariance
            sigma_Q = 1E-13;
            obj.QTilde0 = eye(3, 3) * (sigma_Q^2);
        end
    end
end

