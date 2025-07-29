function [estimObj] = estimation_module(estimObj, OrbitObj, ...
     SAVE, name, filterType)
    %%                        ESTIMATION MODULE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 11/01/2023
    %
    %   Description: This module is in charge of compute the estimation
    %   process using the most suitable filter
    %
    %   Input:
    %       estimObj:  estimation class object
    %       OrbitObj: orbit class object
    %       AttitudeObj: attitude class object
    %       AccObj: acceleration class object
    %       t: time vector
    %       save: save variables option boolean
    %       readData: read data from data files or use current compute data
    %       filterType: type of filter: batch / batch_apriori
    %
    %   Output:
    %       estimObj:  estimation class object
    % --------------------------------------------------------------------%
    % check estimation order
    if(estimObj.EstimData(2) > OrbitObj.OrbitData(9))
        warning("WARNING: trying to estimate a higher harmonic than the " + ...
            "simulated. Not possible");
    elseif(estimObj.EstimData(2) < OrbitObj.OrbitData(9))
        warning("WARNING: trying to estimate a lower harmonic than the " + ...
            "simulated. Not accurate results ");
    end
    disp("Starting estimation process ...");
    
    % ------------------------------------------------------------------- %
    mode = 0;                                         % extra param mode

    % get initial simulation parameters
    OrbitObj = OrbitObj.C_Wmatrix();        % load planet variables
    normalized = OrbitObj.OrbitData(13);    % normalized coefficients
    [estimObj, maxIter] = getParams(estimObj, OrbitObj, normalized, mode);
    
    % load measurements 
    [GG_meas, pos_meas, vel_meas, angVel_meas, ...
    angAcc_meas, B_ACI, TIME] = loadData();

    % include errors
    % % errP = ones(3, length(pos_meas(1, :))) * 0;      % [m] 0.3
    errP = normrnd(0, 0.07, [3, length(pos_meas(1, :))]);
    errV = ones(3, length(vel_meas(1, :))) * 0;    % [m/s] 0.001

    dynamic_meas = [pos_meas + errP;vel_meas + errV;angVel_meas;angAcc_meas];

    % initialize variables
    estimObj = estimObj.initVariables(TIME, maxIter);
    
    % initialize iterator variables
    count = 0;
    error = 100;
    epsilon = estimObj.EstimData(5);                   % convergence value

    tic
    while(abs(error) > epsilon && count < maxIter)
        % filter solver
        if(filterType == 1)
            [XNOT, estimObj] = batchSolver_classic(estimObj, ...
                GG_meas, dynamic_meas, B_ACI, normalized, count);
        elseif(filterType == 2) % CKF
            [estimObj, XNOT] = CKF_solver(estimObj, count, GG_meas, ...
                dynamic_meas,...
                B_ACI, normalized);
        elseif(filterType == 3) % SRIF
            [estimObj, XNOT] = SRIF_solver(estimObj, count, GG_meas, ...
                dynamic_meas,...
                B_ACI, normalized);
        elseif(filterType == 4) % SNC
            [estimObj, XNOT] = SNC_solver(estimObj, count, GG_meas, ...
                dynamic_meas,...
                B_ACI, normalized);
        elseif(filterType == 5) % QR decomposition
            [XNOT, estimObj] = LS_QR_solver(estimObj, ...
                GG_meas, dynamic_meas, B_ACI, normalized, count);
        end
        
        % prefit & postfit per iter
        up = (count + 1)*estimObj.Nm;
        down = up - (estimObj.Nm - 1);
    
        estimObj.prefitj(down:up, :) = estimObj.prefit;
        estimObj.postfitj(down:up, :) = estimObj.postfit;
        
        % compute error.
        error = vecnorm(XNOT);
        disp("XNOT value: " + string(error))
        if(filterType == 1 || filterType == 5)
            error = epsilon;
        end
        
        % update perturbation at epoch.
        estimObj.Xnot = estimObj.Xnot + XNOT;
        estimObj.deltaX = estimObj.deltaX - XNOT;
        
        % update a priori state deviation @ each time step
        estimObj.deltaXj((count + 1)*estimObj.Np-(estimObj.Np-1):(count + 1)*estimObj.Np, :) = ...
            estimObj.pert;
    
        % Increment counter
        count = count + 1;
    end
    toc
    % include last update to the state
    estimObj.Xfilter =  estimObj.Xfilter + XNOT;
    estimObj.Xfilter(1, :) = estimObj.Xfilter(1, :) * estimObj.GM;
    
    % compute error
    [CS_err] = harmonicError(estimObj.Xfilter(1:estimObj.Nc+estimObj.Ns), ...
        OrbitObj.C, OrbitObj.S, ...
        OrbitObj.OrbitData(1), estimObj.Nc, estimObj.Ns);
    disp("SH harmonic error: " + string(vecnorm(CS_err)));
    
    % save data in file 
    if(SAVE == true)
        saveData_estim(estimObj, name, count, CS_err);
    end

% %     % PLOT prefit and postfit over time. XX component
% %     TIME = estimObj.TIME;
% %     figure()
% %     for j =1:count
% %         ind = 6 * j - 5;
% %         subplot(count, 2, j)
% %         plot(TIME, estimObj.prefitj(ind, :), "o", 'Color', 'r')
% %         title('prefit iteration = ' + string(j))
% %     end
% %     i = 1;
% %     for j =count+1:2*count
% %         ind = 6 * i - 5;
% %         subplot(count, 2, j)
% %         plot(TIME, estimObj.postfitj(ind, :), "o", 'Color', 'b')
% %         title('postfit iteration = ' + string(i))
% %         i = i +1;
% %     end
% %     % save compute Cnm and Snm
% %     C_mat = estimObj.Xfilter(1:estimObj.Nc);
% %     S_mat = estimObj.Xfilter(estimObj.Nc+1:estimObj.Nc+estimObj.Ns);
% %     save('C.mat', 'C_mat');
% %     save('S.mat', 'S_mat');
end

