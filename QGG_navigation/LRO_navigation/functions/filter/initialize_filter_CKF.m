function [P0_new, state0_new] = initialize_filter_CKF(time, state0, Y, R0,...
    P0_N, Nt_max, planetParams, Cnm_list, Snm_list, ...
    BN_mat, NB_EARTH, NB_MOON, Q0, Qb, mask)
    
    % state mask (pos, vel & bias)
    mask_state = [ones(6, 1);mask];
    
    % Number of total states
    Nx = length(state0);
    
    % Number of trajectory states (pos & vel)
    Ns = 6;

    PHI0 = reshape(eye(Ns,Ns), [Ns*Ns,1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);

    Pt  = nan(length(time), Nx*Nx);  

    % rotate initial uncertianty to body frame
    P0  = rotate_P0(BN_mat,P0_N);
    BN0 = BN_mat(1:3, :);
    
    Xnot = zeros(Nx, 1); XNOT = Xnot; XNOT_N = XNOT;
    err  = 1; thrs = 1E-15; maxIter = 10; count = 0;
    At = time(2) - time(1);
    while((err > thrs) && (count < maxIter))
        % Integrate trajectory
        X0      = state0(1:Ns) + XNOT_N(1:Ns);
        X0_bias = state0(Ns+1:end) + XNOT(Ns+1:end);
        [~, STATE] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
                 Cnm_list, Snm_list), time(1:Nt_max), [X0;PHI0], options);

        X_N      = STATE(:, 1:Ns)';
        X_bias   = X0_bias.*ones(Ns, Nt_max);
        STM_N    = STATE(:, Ns+1:end);

        % rotate state and STM to body frame
        [X]   = rotate_state(time(1:Nt_max), BN_mat, X_N);
        [STM] = rotate_STM(time(1:Nt_max),BN_mat,STM_N, BN0);

        % process measurements and obtain correction
        for k = 1:Nt_max
            if(k == 1)
                X_hat = -XNOT;
                P = P0;
            end

            % compute previous STM
            if( k > 1)
                PHI_p = reshape(STM(k - 1, :), [Ns, Ns]);
                
                % STM @ ti 2 ti-1
                PHI_ij = reshape(STM(k, :), [Ns,Ns])/PHI_p;
            else
                % STM @ ti 2 ti-1
                PHI_ij = reshape(STM(k, :), [Ns,Ns]);
            end

            % compute measurements and partials
            maxInd         = 3 * k; minInd = maxInd - 2;
            BN             = BN_mat(minInd:maxInd, :);
            BODYMOON_J2000 = NB_MOON(minInd:maxInd, :)';
            
            X_N_k          = blkdiag(BN', BN') * X(:, k);
            [Yc, Hi, ~]    = compute_measurements_filter(planetParams, ...
                        time(k), X_N_k', Cnm_list, Snm_list, BN, BODYMOON_J2000);
            dY             = Y(:, k) - Yc - X_bias(:, k);
            
            % Include process noise
            Q_N = processNoise(Q0, Qb, At, Nx);
            Q  = rotate_processNoise(BN_mat,Q_N, k);

            % Augment STM to include the bias
            PHI_tot = [PHI_ij, zeros(6,6);...
                       zeros(6,6), eye(6)];

            % apply measurement mask
            dY_used = dY(logical(mask));
            Hi_used = Hi(logical(mask), logical(mask_state));
            R0_used = R0(logical(mask), logical(mask));
            
            % apply state mask
            PHI_used = PHI_tot(logical(mask_state), logical(mask_state));
            Q_used   = Q(logical(mask_state), logical(mask_state));
            P_used   = P(logical(mask_state), logical(mask_state));
            X_hat_used    = X_hat(logical(mask_state));

            % run CKF
            [X_hat_new, P_new] = CKF(dY_used, Hi_used, R0_used, P_used, ...
                X_hat_used, PHI_used, Q_used);

            % update states inside the mask
            X_hat(logical(mask_state))                      = X_hat_new;
            P(logical(mask_state), logical(mask_state))     = P_new;

            % update current state
            X(:, k)      = X(:, k)      + X_hat(1:Ns);
            X_bias(:, k) = X_bias(:, k) + X_hat(Ns+1:end);

            % % err = abs(X(:, k) - state_true(k, 1:6)');
    
            % current uncertainty
            Pt(k, :) = reshape(P, [1, Nx*Nx]);
        end
        
        % restart filter
        PHI          = reshape(STM(end, :), [Ns, Ns]);
        PHI_end_init = [PHI, zeros(6,6);...
                       zeros(6,6), eye(6)];
        Xnot         = PHI_end_init\X_hat;

        % add corrections (body frame)
        XNOT   = XNOT + Xnot;
        BN0    = BN_mat(1:3, :);
        XNOT_N = blkdiag(BN0', BN0', eye(6)) * XNOT; 

        err    = vecnorm(Xnot);
        count  = count + 1;
    end
    % state in body frame
    Xf_B = [X(:, 1);X_bias(:, 1)];
    Pf_B = reshape(Pt(1, :), [Nx, Nx]);
    
    % return state in inertial frame
    BN0  = BN_mat(1:3, :);
    A    = blkdiag(BN0', BN0', eye(6)); 
    P0_new      = A * Pf_B * A.';
    state0_new  = A * Xf_B;
end

