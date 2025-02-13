function [varObj, XNOT] = SNC_solver(varObj, count, GG_meas, dynamic_meas,...
    B_ACI_mat, normalized)
    % SH maximum degree
    n_max  = varObj.n_max;

    % get dynamical measurements
    ri = dynamic_meas(1:3, :);
    %vi = dynamic_meas(4:6, :);
    angVel = dynamic_meas(7:9, :);
    angAcc = dynamic_meas(10:12, :);
    
    % get pole parameters. [W, W0, RA, DEC]
    W = varObj.poleParams(1);
    W0 = varObj.poleParams(2);
    RA = varObj.poleParams(3);
    DEC = varObj.poleParams(4);

    % Initial state @ current iteration
    X = varObj.X0 + varObj.Xnot;

    % Compute state dynamics. ODE 113
    STM = reshape(eye(varObj.Np, varObj.Np), [varObj.Np^2, 1]);
    ref_vec = [X; STM];

    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [~, state] = ode113(@(t, x) dynamic_propagator(t, x, varObj), ...
        varObj.TIME, ref_vec, options);
    
    varObj.Xfilter = state(:, 1:varObj.Np)';
    varObj.Xfilterj((count + 1)*varObj.Np-(varObj.Np-1):(count + 1)*varObj.Np, :) = ...
        varObj.Xfilter;
    
    STM = state(:, varObj.Np+1:end)';

    % get SH matrices at nominal
    [C_mat, S_mat] = get_SHmatrices(varObj, X);
    GM = varObj.GM * X(1);

    % go at each time
    for k =1:varObj.Nt
        % initial CKF and EKF estimations
        if(k == 1)
            X_hat = varObj.deltaX;
            P = varObj.P0;
        end

        % rotation matrices @ current time
        up = 3*k;
        down = up -2;
        B_ACI = B_ACI_mat(down:up, :);
        
        % ACI 2 ACAF rotation matrix
        Wt = W0 + W * varObj.TIME(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        B_ACAF = B_ACI * ACAF_ACI';
    
        % state vector @ current time. ACAF
        r_ACAF = ACAF_ACI * ri(:, k);
        %v_ACAF = ACAF_ACI * vi(:, k);
        
        % state vector @ current time. ACI
        %r_ACI = ri(:, k);
        %v_ACI = vi(:, k);

        % current angular vel / acc
        omega_k = angVel(:, k);
        Omega_k = angAcc(:, k);

        % compute cross product operator
        [omega2M, OmegaM] = angularMatrix(omega_k, Omega_k);

        % STM @ current time. ti 2 t0
        PHI = reshape(STM(:, k), [varObj.Np, varObj.Np]);
        PHI = PHI(1:varObj.Np, 1:varObj.Np);
    
        % STM @ previous time. ti-1 2 t0
        if( k > 1)
            PHI_p = reshape(STM(:, k - 1), [varObj.Np, varObj.Np]);
            PHI_p = PHI_p(1:varObj.Np, 1:varObj.Np);
            
            % STM @ ti 2 ti-1
            PHI_ij = PHI/PHI_p;
        else
            % STM @ ti 2 ti-1
            PHI_ij = PHI;
        end

          % compute visibility matrix
         hi_9t = potentialGradient_Cnm(n_max, r_ACAF, ...
            varObj.Re, GM, B_ACAF, normalized);
         hg_t = [hi_9t(1, :); hi_9t(4, :); hi_9t(7, :); ...
             hi_9t(5, :); hi_9t(8, :); hi_9t(9, :)];

          % compute visibility matrix. Extra paramters
         [he_t, b] = H_extraParam(varObj, k);

         % add both
         if(varObj.extraParam~=0)
             Hi_t = [hg_t, he_t];
         else
             Hi_t = hg_t;
         end

         % compute observable & meas @ current time
         [~, ~, ddU_ACAF] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                       r_ACAF, varObj.Re, GM, ...
                                       normalized);
         ddU_B = B_ACAF * ddU_ACAF * B_ACAF';
        
         yc  = ddU_B - omega2M - OmegaM + b;
         Yc = [yc(1,1); yc(1, 2); yc(1,3); ...
             yc(2,2); yc(2,3); yc(3,3)];   
         
         yi = -GG_meas(down:up, :);
         Yi = [yi(1,1); yi(1,2); yi(1,3);...
             yi(2,2); yi(2,3); yi(3,3)];

         dY = Yi - Yc;

         % reshape
         dY = reshape(dY, [varObj.Nm, 1]);
        
         % update process noise
         if(k == 1)
            Q = zeros(varObj.Np, varObj.Np);
        else
            At = varObj.TIME(k) - varObj.TIME(k - 1);
            QT = varObj.QTilde0;
            Q = processNoise(QT, At, varObj);
        end

        % run filter
        [X_hat, P, ~] = CKF_SNC(dY, Hi_t, varObj.R0, P, X_hat, PHI_ij, Q);

        % covariance matrix
        varObj.Pj(k*varObj.Np-(varObj.Np - 1):k*varObj.Np, :, count+1) = P;
        varObj.sigma_t(:, k) = sqrt(diag(P));

        % state pertubation estimation
        varObj.pert(:, k) = X_hat;
            
        % compute prefit @ current time
        varObj.prefit(:, k) = dY;
        
        % compute postfit @ current time
        varObj.postfit(:, k) = dY - (Hi_t * X_hat);
    end 

    % map to initial state
    PHI = reshape(STM(:, varObj.Nt), [varObj.Np, varObj.Np]);

    % compute XNOT.
    XNOT = PHI\X_hat;

    % compute a posteriori covariance
    varObj.Pj0(:, :, count+1) = inv(PHI) * P * inv(PHI)';
end