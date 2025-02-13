function [XNOT, P, estimObj] = batchSolver_apriori(estimObj, OrbitObj, ...
                                    t, P0, Xnot_j)
    %%                      BATCH SOLVER APRIORI FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 05/02/2023
    %
    %   Description: Funtion with the batch routine. 
    %
    %   Input: estimObj: estimation object
    %          t: time vector
    %          P0: a priori covariance
    %          Xnot_j: initial a prioiri deviation, iteration j
    %
    %   Output: XNOT: state deviation after processing meas
    %           P: covariance matrix, state vector
    %           estimObj: estimation data object
    %
    % --------------------------------------------------------------------%

    % get planet pole parameters
    W = OrbitObj.W;
    W0 = OrbitObj.W0;
    RA = OrbitObj.RA;
    DEC = OrbitObj.DEC;

    % get normalized harmonics config
    Normalized = OrbitObj.OrbitData(13);

    % Weights
    sigma2_ii = estimObj.EstimData(3)^2;
    sigma2_ij = estimObj.EstimData(4)^2;
    Ri = diag([sigma2_ii, sigma2_ij, sigma2_ij, sigma2_ij, sigma2_ii, ...
        sigma2_ij, sigma2_ij, sigma2_ij, sigma2_ii]);

    % Init normal equation values
    Ax = inv(P0);
    Nx = P0\(-Xnot_j);

    % time loop
    for k = 1:estimObj.Nt
        % ACI 2 ACAF rotation matrix
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % current position. ACI frame
        r_ACI = estimObj.ri(:, k);
        r_ACAF = ACAF_ACI * r_ACI;

        % current angular vel / acc
        omega_k = estimObj.angVel(:, k);
        Omega_k = estimObj.angAcc(:, k);

        % compute cross product operator
        [omega2M, OmegaM] = angularMatrix(omega_k, Omega_k);

        % rotation matrix ACAF 2 body @ current time
        up = 3*k;
        down = up - 2;
        B_ACI = estimObj.BN(down:up, :);
        B_ACAF = B_ACI * ACAF_ACI';

        % construct measurements tensor. Body frame
        up = estimObj.Nm^0.5 * k;
        down = up - (estimObj.Nm^0.5 - 1);
        
        accT = estimObj.accT(down:up, :);
        Odata = reshape(-accT, [estimObj.Nm, 1]);
        
        % compute tensor at nominal
        [~, ~, ddU_ACAF] = potentialGradient_nm(estimObj.C, estimObj.S, ...
            estimObj.EstimData(2), r_ACAF, estimObj.Re, estimObj.X(1), ...
            Normalized);
        ddU_B = B_ACAF * ddU_ACAF * B_ACAF';

        % compute tensor at nominal
        Cdata = ddU_B - omega2M - OmegaM;
        Cdata = reshape(Cdata, [estimObj.Nm, 1]);

        % O - C value
        deltaY = Odata - Cdata;
        
        % compute vibility matrix. At nominal
        H = potentialGradient_Cnm(estimObj.EstimData(2), r_ACAF, ...
            estimObj.Re, estimObj.X(1), B_ACAF, Normalized);

        % Batch filter
        [Ax, Nx] = batchFilter(deltaY, H, Ri, Ax, Nx);

        % prefit 
        estimObj.pref(:, k) = deltaY;

        % store H matrix
        up = estimObj.Nm * k;
        down = up - (estimObj.Nm - 1);
        estimObj.H(down:up, :) = H;

        % variance over time
        estimObj.Pt(:, k) = diag(inv(Ax));
    end
    
    % solve Normal equation
    XNOT = Ax\Nx;

    % compute covariance estimation
    P = inv(Ax);
end

