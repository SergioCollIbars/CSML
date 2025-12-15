function [X, Px, Xc, P_grav, h, ax] = do_SLAM_v3(metaData_file, time, ...
    state_true, planetParams,poleParams, instrumentParams, ...
    BN_mat, Cnm_t, Snm_t, Y)
    
    % load filter parameters
    [R0, P0,  P0c_grav, Q0, Qb, delta_state0, Cnm, Snm, Nbatch] = ...
        loadFilterParams(metaData_file, planetParams, instrumentParams,...
        Cnm_t, Snm_t);

    state0 = delta_state0+state_true(1, 1:12)';

    % state, time numbers
    [Ncoeff, Nsoeff, Ncs]  = count_num_coeff(planetParams(3)); 
    n_solve = planetParams(3);
    Ns = 12; Nt = length(time);
    At = time(2) - time(1);
    Nx = 12 + Ncs;

    % initiate filter
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    X = zeros(Ns, Nt);
    X(:, 1) = state0; X0 = X(:, 1);
    
    Pxc = zeros(Ns, Ncs); Pcx = zeros(Ncs, Ns);
    Pcc = diag(diag(P0c_grav));
    P = [P0, Pxc;...
        Pcx, Pcc];

    Xc_t   = mat2list(Cnm_t, Snm_t, Ncoeff, Nsoeff);
    Xc_0   = mat2list(Cnm, Snm, Ncoeff, Nsoeff);
    Xc     = Xc_0; 

    % covariance & state correction at each time
    Px   = ones(Nt, Ns*Ns) * NaN;
    Px(1, :) = reshape(P0, [Ns*Ns, 1]);
    Pc_grav  = P0c_grav; 

    errNorm_pos   = nan(1, Nt); errNorm_vel   = nan(1, Nt);
    errNorm_bias  = nan(6, Nt);
    std_pos       = nan(1, Nt); std_vel       = nan(1, Nt);
    std_bias      = nan(6, Nt);

    % start real-time plots
    [fig, h, ax] = start_plots(time, planetParams(3), Xc_t, ...
        Xc_0, P0c_grav, Nbatch);
    pause(1)

    % print once (on a new line, doesn't erase anything)
    fprintf('\nProgress:   0%%'); count = 0; flag_count = 0;
    for k = 2:Nt
        % initial conditions
        t_span = [time(k - 1), time(k)];
        
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
        [X_hat, P] = EKF(dY, [Hi, Hc], R0, P, PHI_tot, Q, Ns, Ncs);

% %         ee = sum(eig(P) < 0);
% %         if(ee > 0)
% %             disp(eig(P));
% %         end

        % update nominal
        X(:, k) = state' + X_hat(1:Ns);
        X0 = X(:, k);

        % save covariance
        Px(k, :) = reshape(P(1:Ns, 1:Ns), [1, Ns*Ns]);

        if mod(k,Nbatch) == 0
           flag = 1;     % trigger
           flag_count = flag_count + 1;
        else
           flag = 0;     % reset
           count = count + 1;
        end

        % update the percentage in-place (4 characters: '100%')
        pct = 100 * count / Nbatch;
        fprintf('\b\b\b\b%3.0f%%', pct);
        if(flag)
            % calibrated signal
            if(flag_count == 1)
                st = 10;
            else
                st = Nbatch - 50;
            end
            % % Y_cal   = Y(:, st:k) - X(7:end, st:k);
            Y_cal   = Y(:, st:k);

            [Xg_smooth, Pg_smooth, NSM_err] = gravEstim_smoothing(time(st:k), Y_cal, ...
                planetParams, poleParams, X(:, st:k), P0c_grav, ...
                Cnm, Snm, R0, Px(st:k, :));
            fprintf('\n   Gravity estimation number %d. Initial error: %.2e. Final error: %.2e.\n', ...
            flag_count, NSM_err(1), NSM_err(end));
            fprintf('   Running smoothing EKF ...')

            % update gravity field
            [Cnm_new, Snm_new] = list2mat(n_solve, Ncoeff, Nsoeff, Xg_smooth);

            Cnm(1:n_solve+1,1:n_solve+1) = Cnm_new;
            Snm(1:n_solve+1,1:n_solve+1) = Snm_new;

            % assign initial values (start at final time)
            P0_smooth              = reshape(Px(1, :), [Ns, Ns]);
            Pc_grav(2:end, 2:end)  = Pg_smooth;
    
            % run EKF
            [Xs_smooth, Ps_smooth] = EKF_smoothing(time(1:k), planetParams, ...
                poleParams, Y(:, 1:k), X(:, 1), P0_smooth, ...
                Pc_grav, BN_mat, Cnm, Snm, Q0, Qb, R0);
            
            % update parameters
            p = reshape(Ps_smooth(end, :), [Nx, Nx]);
            % % P = [p,Pxc;Pcx, Pc_grav];
            P = p;
            for j = 1:k
                b = reshape(Ps_smooth(j, :), [Nx, Nx]);
                Px(j, :) = reshape(b(1:Ns, 1:Ns), [Ns*Ns, 1]);
            end
            X(:, 1:k)  = Xs_smooth;
            X0         = X(:, k);
            
            Xc         = Xg_smooth;

            % compute gravity coefficient error
            errC   = abs(Xc - Xc_t);
            sigmaC = diag(sqrt(Pc_grav));

            % Update plots
            [errNorm_pos, errNorm_vel, ...
            errNorm_bias, std_pos, std_vel, std_bias] = ...
            update_plots(time, k, state_true, X, Px, h, ax, errNorm_pos,...
            errNorm_vel, errNorm_bias, std_pos, std_vel, std_bias, ...
            errC, sigmaC, planetParams(3));
            
            fprintf('\nDone.\n');
            fprintf('\nProgress:   0%%'); count = 0;
        end
    
        % stop if figure is closed
        if ~ishandle(fig)
            disp('Figure closed, stopping simulation.');
            break;
        end
    end

    % save final gravity uncertainty
    P_grav = zeros(Ncs, Ncs);
    P_grav(2:end, 2:end) = Pg_smooth;
end