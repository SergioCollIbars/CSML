function [P0_new, state0_new] = initialize_filter_CKF(time, state0, Y, R0,...
    P0, Pc, Nt_max, planetParams, poleParams, Cnm, Snm, BN_mat, Q0, Qb)
    
    Ns = length(state0);
    [~, ~, Ncs]  = count_num_coeff(planetParams(3)); 
    Nxx = Ncs + Ns;

    PHI0 = reshape(eye(Nxx,Nxx), [Nxx*Nxx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % consider parameters
    Pcc = Pc;
    Pxc = zeros(Ns, Ncs);
    Pcx = zeros(Ncs, Ns);

    P0  = [P0, Pxc; Pcx, Pcc]; 
    Pt  = nan(length(time), Nxx*Nxx);  

    Xnot = zeros(Nxx, 1); XNOT = Xnot;
    err  = 1; thrs = 1E-15; maxIter = 10; count = 0;
    At = time(2) - time(1);
    while((err > thrs) && (count < maxIter))
        % run Batch filter
        X0 = state0 + XNOT(1:Ns);
        [~, STATE] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
            Cnm, Snm), time(1:Nt_max), [X0;PHI0], options);

        X    = STATE(:, 1:Ns)';
        STM  = STATE(:, Ns+1:end);

        for k = 1:Nt_max
            if(k == 1)
                X_hat = -XNOT;
                P = P0;
            end

            % compute previous STM
            if( k > 1)
                PHI_p = reshape(STM(k - 1, :), [Nxx, Nxx]);
                
                % STM @ ti 2 ti-1
                PHI_ij = reshape(STM(k, :), [Nxx,Nxx])/PHI_p;
            else
                % STM @ ti 2 ti-1
                PHI_ij = reshape(STM(k, :), [Nxx,Nxx]);
            end


            % compute measurements and partials
            maxInd      = 3 * k; minInd = maxInd - 2;
            BN          = BN_mat(minInd:maxInd, :);
            [Yc, Hi, Hc, ~]    = compute_measurements_filter(planetParams, ...
                poleParams, time(k), X(:, k)', Cnm, Snm, BN);
            dY          = Y(:, k) - Yc;
            
            % Include process noise
            Q = processNoise(Q0, Qb, At, Nxx);

            % run CKF
            [X_hat, P] = CKF(dY, ...
                 [Hi, Hc], R0, P, X_hat, PHI_ij, Q, Ns, Ncs);

            % update current state
            X(:, k) = X(:, k) + X_hat(1:Ns);
    
            % current uncertainty
            Pt(k, :) = reshape(P, [1, Nxx*Nxx]);
        end
        
        % restart filter
        PHI_end_init = reshape(STM(end, :), [Nxx, Nxx]);
        Xnot = PHI_end_init\X_hat;

        % add corrections
        XNOT = XNOT + Xnot;

        err    = vecnorm(Xnot);
        count  = count + 1;
    end

    P0_new      = reshape(Pt(1, :), [Nxx, Nxx]);
    state0_new  = X(:, 1);
end

