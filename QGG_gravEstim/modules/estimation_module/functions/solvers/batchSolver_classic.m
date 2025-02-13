function [XNOT, estimObj] = batchSolver_classic(estimObj,...
    GG_meas, dynamic_meas, B_ACI_mat, normalized , count)
    %%                      BATCH SOLVER CLASSIC FUNCTION
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

    % Weights
    Ri = estimObj.R0;

    % Init normal equation values
    Ax = 0;
    Nx = 0;
    
    % time vector
    t = estimObj.TIME;

    % GM value
    GM = estimObj.GM;
    
    % time step
    At = t(2) - t(1);

    % visibility matrix at each time step
    Ht = ones(estimObj.Nm * estimObj.Nt, estimObj.Np)*NaN;
    % Obtained
    B = [-1.41000139962798e-07;...
        4.42692393515439e-09;...
        8.21383297049534e-09;...
        -3.18001688021536e-07;...
        -9.74999439415953e-06;...
        -5.65000358162869e-08];

    D = [ 3.78130353377026e-17;...
        8.11853407415446e-17;...
        -6.23503707619306e-18;...
        1.13081945673253e-15;...
        5.41557860429089e-15;...
        1.96147289639821e-17];

% %     %        Required
% %     B = [-1.40999980249042e-07;...
% %         4.44012601351793e-09;...
% %         8.20964327520085e-09;...
% %         -3.18000392726029e-07;...
% %         -9.75000015292926e-06;...
% %         -5.64996665239974e-08];
% % 
% %     D = [3.75113030820805e-17;...
% %         5.94765447757001e-17;...
% %         7.04643715521234e-19;...
% %         1.12864253647863e-15;...
% %         5.42507893184713e-15;...
% %         1.90770566034683e-17];
      
% %        B = zeros(6, 1);
% %        D = B;

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

        deltaY = Ydata + B + k * D * At;
        
        % compute vibility matrix for the SH coefficients
        hi_9t = potentialGradient_Cnm(estimObj.EstimData(2), r_ACAF, ...
            estimObj.Re, GM, B_ACAF, normalized);
        hg_t = [hi_9t(1, :); hi_9t(4, :); hi_9t(7, :); ...
             hi_9t(5, :); hi_9t(8, :); hi_9t(9, :)];

        % compute visibility matrix. Extra parameters
        [he_t, ~, ~] = H_extraParam(estimObj, k, At);

         % add both
         if(estimObj.extraParam~=0)
             H = [hg_t, he_t];
         else
             H = hg_t;
         end

        % Batch filter
        [Ax, Nx] = batchFilter(deltaY, H, Ri, Ax, Nx);

        % prefit 
        estimObj.prefit(:, k) = deltaY;

        % save H matrix
        ind1 = k * estimObj.Nm;
        ind2 = ind1 - (estimObj.Nm - 1);
        Ht(ind2:ind1, :) = H;
    end

    % solve Normal equation
    XNOT = Ax\Nx;
    
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

