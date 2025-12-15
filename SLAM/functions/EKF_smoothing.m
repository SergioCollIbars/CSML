function [X_smooth, P_smooth] = EKF_smoothing(time, planetParams, poleParams, Y, X0, P0, ...
    P0c_grav, BN_mat, Cnm, Snm, Q0, Qb, R0)
    
    % Initialize filter
    Nstart = 10;
    [P0_new, state0_new] = initialize_filter_CKF(time(1:Nstart), X0, Y, R0,...
        P0, P0c_grav, Nstart, planetParams, poleParams, Cnm, Snm, BN_mat, Q0, Qb);

    X0 = state0_new;
    
    Ns = length(X0);
    [~, ~, Ncs]  = count_num_coeff(planetParams(3)); 
    Nx = Ncs + Ns;

    % initiate filter
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % output variables
    X_smooth = ones(Ns, length(time)) * NaN;
    P_smooth = nan(length(time), Nx*Nx);
    
    % assign initial values
    P_old = P0_new;
    P_smooth(1, :) = reshape(P_old, [Nx*Nx, 1]);
    X_smooth(:, 1) = X0; 
    
    At = time(2) - time(1);
    fprintf('  0%%');
    for k = 2:length(time)
        pct = 100 * k / length(time);
        fprintf('\b\b\b\b%3.0f%%', pct);

        % initial conditions
        t_span = [time(k-1), time(k)];
        
        % Compute dynamics. ODE 113
        [~, STATE] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
            Cnm, Snm), t_span, [X0;PHI0], options);
        
        state   = STATE(end, 1:Ns);
        PHI_tot = reshape(STATE(end, 13:end), [Nx, Nx]);
        
        % compute predicted measurements (TBD)
        maxInd      = 3 * k; minInd = maxInd - 2;
        BN          = BN_mat(minInd:maxInd, :);
        [Yc, Hi, Hc, ~]    = compute_measurements_filter(planetParams, ...
            poleParams, time(k), state, Cnm, Snm, BN);
        dY          = Y(:, k) - Yc; 
    
        % Include process noise
        Q = processNoise(Q0, Qb, At, Nx);
    
        % Run EKF process
        [X_hat, P_new] = EKF(dY, [Hi, Hc], R0, P_old, PHI_tot, Q, Ns, Ncs);
    
    % %             ee = sum(eig(P) < 0);
    % %             if(ee > 0)
    % %                 disp(eig(P));
    % %             end
    
        % update nominal & save states
        X_smooth(:, k) = state' + X_hat(1:Ns);
        X0 = X_smooth(:, k);   
    
        P_smooth(k, :) = reshape(P_new, [Nx*Nx, 1]);
        P_old          = P_new;
    end
    fprintf('\n');
end

