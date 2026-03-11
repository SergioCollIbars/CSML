function [X, Pt, Xhat, Xnot, pref, posf] = CKF_solver_EPHEM(TIME, X0, ...
    Xnot, P0, R0, Q0, Qb, meas, planetParams, BN_matrix, C_mat, ...
    S_mat,posE, posM, posS, measMask)
    % Run CKF solver using only the EPHEMERIDES dynamics.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of parameters
    Nt = length(TIME);
    Ns = 12;

    % state maks
    stateMask = [ones(6, 1);measMask];

    % state values
    X0 = X0 + Xnot;
    STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
    initState = [X0; STM0];

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, STATE] = ode113(@(t, x) EOM_navigation_EPHEM(t, x, planetParams, ...
        C_mat, S_mat), TIME, initState, options);
    X    = STATE(:, 1:Ns)';
    STM  = STATE(:, Ns+1:end);

    Pt = ones(Nt, Ns*Ns) * NaN;
    Xhat = ones(Ns, Nt) * NaN;

    At = (TIME(2) - TIME(1)); % [-]
    
    % start filter
    for j = 1:Nt
        % initial CKF and EKF estimations
        if(j == 1)
            X_hat = -Xnot;
            P = P0;
        end
        
        % compute previous STM
        if( j > 1)
            PHI_p = reshape(STM(j - 1, :), [Ns, Ns]);
           
            % STM @ ti 2 ti-1
            PHI_ij = reshape(STM(j, :), [Ns,Ns])/PHI_p;
        else
            % STM @ ti 2 ti-1
            PHI_ij = reshape(STM(j, :), [Ns,Ns]);
        end

         % compute gravity tensor measurements, accelerometer measurements
         % and partials
        state = X(:, j)';
        maxInd = 3 * j; minInd = maxInd -2;
        BN = BN_matrix(minInd:maxInd, :);
        [T_c, ~] = compute_measurement_EPHEM(TIME(j), state, planetParams,...
            BN, C_mat, S_mat, 0, posE(:, j), posM(:, j), posS(:, j));
        [Hmeas] = compute_meas_partials_EPHEM(TIME(j), state, planetParams,...
                    BN, C_mat, S_mat, posE(:, j), posM(:, j), posS(:, j));
        
        % measurement residuals
        dY = meas(:, j) - T_c;
        
        % process noise
        Q = processNoise(Q0, 0, At, Qb, "SNC", Ns);
        if(j== 1)
            Q = Q.*0;
        end
        
        % apply mask
        dY_used    = dY(logical(measMask));
        Hmeas_used = Hmeas(logical(measMask), logical(stateMask));
        R0_used    = R0(logical(measMask), logical(measMask));
        Q_used     = Q(logical(stateMask), logical(stateMask));
        PHI_used   = PHI_ij(logical(stateMask), logical(stateMask));
        P_used     = P(logical(stateMask), logical(stateMask));
        X_hat_used = X_hat(logical(stateMask));

        % run CKF
        [X_hat_new, P_new, ~, ~] = CKF(dY_used, ...
             Hmeas_used, R0_used, P_used, X_hat_used, PHI_used, Q_used);

        Xhat(logical(stateMask), j) = X_hat_new;
        X_hat(logical(stateMask))   = X_hat_new;

        % update current state
        X(logical(stateMask), j) = X(logical(stateMask), j) + X_hat_new;

        % current uncertainty
        P_sym = (P_new + P_new.')/2;
        P(logical(stateMask), logical(stateMask)) = P_sym;
        Pt(j, :) = reshape(P, [1, Ns*Ns]);
    end
    
    % restart filter
    PHI0 = reshape(STM(end, :), [Ns, Ns]);
    Xnot = PHI0\X_hat;

    % compute prefit and postfit
    Nm = length(R0);
    pref = ones(Nm, Nt) * NaN;
    posf = ones(Nm, Nt) * NaN;
    for j = 1:Nt
        PHI_i0 = reshape(STM(j, :), [Ns,Ns]);

        state = X(1:Ns, j)';
        maxInd = 3 * j; minInd = maxInd -2;
        BN = BN_matrix(minInd:maxInd, :);
        [T_c, ~] = compute_measurement_EPHEM(TIME(j), state, planetParams,...
            BN, C_mat, S_mat, 0, posE(:, j), posM(:, j), posS(:, j));
        [Hi] = compute_meas_partials_EPHEM(TIME(j), state, planetParams,...
                    BN, C_mat, S_mat, posE(:, j), posM(:, j), posS(:, j));

        pref(:, j) = meas(:, j) - T_c;
        posf(:, j) = pref(:, j) - Hi * (PHI_i0 * Xnot);
    end
end