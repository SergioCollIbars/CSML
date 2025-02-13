function [XNOT, estimObj] = LS_QR_solver(estimObj,...
    GG_meas, dynamic_meas, B_ACI_mat, normalized , count)
    %%                     LS QR SOLVER FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 05/02/2023
    %
    %   Description: Funtion with the batch routine. 
    %
    %   Input: estimObj: estimation object
    %          estimObj: estimation data object
    %          OrbitObj: orbit data object 
    %
    %   Output: XNOT: state estimation vector
    %           estimObj: estimation data object
    %
    % --------------------------------------------------------------------%
    % get planet pole parameters
    W = estimObj.poleParams(1);
    W0 = estimObj.poleParams(2);
    RA = estimObj.poleParams(3);
    DEC = estimObj.poleParams(4);

    % get dynamical measurements
    ri = dynamic_meas(1:3, :);
    %vi = dynamic_meas(4:6, :);
    angVel = dynamic_meas(7:9, :);
    angAcc = dynamic_meas(10:12, :);

    % time vector
    t = estimObj.TIME;

    % Weights
    Wmeas = inv(estimObj.R0);

    % GM value
    GM = estimObj.GM;
    
    % time step
    At = t(2) - t(1);

    % visibility matrix at each time step
    Ht = ones(estimObj.Nm * estimObj.Nt, estimObj.Np)*NaN;
    
    % measurement matrix at each time step
    Yt = ones(estimObj.Nm * estimObj.Nt, 1)*NaN;

    % Information matrix. Ax 
    Ax  = 0;

     % time loop
    for k = 1:estimObj.Nt
        % ACI 2 ACAF rotation matrix
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % current position. ECI frame
        r_ACI = ri(:, k);
        r_ACAF = ACAF_ACI * r_ACI;

        % current angular vel / acc
        omega_k = angVel(:, k);
        Omega_k = angAcc(:, k);

        % compute cross product operator
        [omega2M, OmegaM] = angularMatrix(omega_k, Omega_k);

        % rotation matrix ECEF 2 body @ current time
        up = 3*k;
        down = up - 2;
        B_ACI = B_ACI_mat(down:up, :);
        B_ACAF = B_ACI * ACAF_ACI';

        % construct measurements tensor
        accT = -GG_meas(down:up, :) + omega2M + OmegaM;
        
        % get data measurements
        Ydata = [accT(1,1); accT(1, 2); accT(1, 3); accT(2, 2); ...
            accT(2, 3); accT(3,3)];
        deltaY = Ydata;
        
        % compute vibility matrix
        hi_9t = potentialGradient_Cnm(estimObj.EstimData(2), r_ACAF, ...
            estimObj.Re, GM, B_ACAF, normalized);
        hg_t = [hi_9t(1, :); hi_9t(4, :); hi_9t(7, :); ...
             hi_9t(5, :); hi_9t(8, :); hi_9t(9, :)];

        % compute visibility matrix. Extra paramters
        [he_t, ~, ~] = H_extraParam(estimObj, k, At);

         % add both
         if(estimObj.extraParam~=0)
             H = [hg_t, he_t];
         else
             H = hg_t;
         end

        % prefit 
        estimObj.prefit(:, k) = deltaY;
        
        % save Y matrix
        ind1 = k * estimObj.Nm;
        ind2 = ind1 - (estimObj.Nm - 1);
        Yt(ind2:ind1, :) = deltaY;

        % save H matrix
        Ht(ind2:ind1, :) = H;

        % Ax
        Ax = Ax + H' * Wmeas * H;
    end

    % solve QR decomposition
    XNOT = fixed.qrMatrixSolve(Ht,Yt);
    
    % compute postfit
    for k = 1:estimObj.Nt
        ind1 = k * estimObj.Nm;
        ind2 = ind1 - (estimObj.Nm - 1);
        estimObj.postfit(:, k) = estimObj.prefit(:, k) - ...
            Ht(ind2:ind1, :) * XNOT;
    end

    % compute covariance estimation
    estimObj.Pj0(:, :, count+1) = inv(Ax);
end

