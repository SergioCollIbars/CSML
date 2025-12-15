function [X, Px, Xc, P0c_grav, h, ax] = do_SLAM(time, state_true, planetParams, ...
    poleParams, instrumentParams, BN_mat, Cnm_t, Snm_t, Y)
    
    % load filter parameters
    [R0, P0,  P0c_grav, Q0, Qb, delta_state0, Cnm, Snm, Nbatch, ~] = ...
        loadFilterParams(planetParams, instrumentParams, Cnm_t, Snm_t);

    state0 = delta_state0+state_true(1, 1:12)';
    state0(7:end) = zeros(6, 1);

    % state, time numbers
    [Ncoeff, Nsoeff, Ncs]  = count_num_coeff(planetParams(3)); 
    Ns = 12; Nt = length(time); Nm = length(R0);
    At = time(2) - time(1);
    Nx = 12 + Ncs;

    % initiate filter
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    X = zeros(Ns, Nt);
    X(:, 1) = state0; X0 = X(:, 1);
    
    Pxc = zeros(Ns, Ncs); Pcx = zeros(Ncs, Ns);
    Pcc = diag(diag(P0c_grav)); P0c = diag(diag(P0c_grav));
    P = [P0, Pxc;...
        Pcx, Pcc];

    Xc_t   = mat2list(Cnm_t, Snm_t, Ncoeff, Nsoeff);
    Xc_0   = mat2list(Cnm, Snm, Ncoeff, Nsoeff);
    Xc     = Xc_0; 

    % covariance & state correction at each time
    Px   = ones(Nt, Ns*Ns) * NaN;
    Xhat = ones(Ns, Nt) * NaN;
    Px(1, :) = reshape(P0, [Ns*Ns, 1]);
    
    % prefit and postfit
    pref = ones(Nm, Nt) * NaN;
    posf = ones(Nm, Nt) * NaN;

    errNorm_pos   = nan(1, Nt); errNorm_vel   = nan(1, Nt);
    errNorm_bias  = nan(6, Nt);
    std_pos       = nan(1, Nt); std_vel       = nan(1, Nt);
    std_bias      = nan(6, Nt);

    % start real-time plots
    [fig, h, ax] = start_plots(time, planetParams(3), Xc_t, ...
    Xc_0, P0c_grav, Nbatch);
    
    total_inf = inv(P0c_grav(2:end, 2:end)); kmin = 1;
    flag_count= 0; 
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
        [Yc, Hi, Hc, Hp]    = compute_measurements_filter(planetParams, ...
            poleParams, time(k), state, Cnm, Snm, BN);
        dY          = Y(:, k) - Yc; 
        pref(:, k)  = dY;

        % Include process noise
        Q = processNoise(Q0, Qb, At, Nx);

        % Run EKF process
        [X_hat, P] = EKF(dY, [Hi, Hc], R0, P, PHI_tot, Q, Ns, Ncs);

        ee = sum(eig(P) < 0);
        if(ee > 0)
            disp(eig(P));
        end

        % update nominal
        X(:, k) = state' + X_hat(1:Ns);
        X0 = X(:, k);
        
        % postfit
        posf(:, k) = pref(:, k) - Hi * X_hat(1:Ns);

        % save covariance
        Px(k, :) = reshape(P(1:Ns, 1:Ns), [1, Ns*Ns]);

        % save correction
        Xhat(:, k) = X_hat(1:Ns);

        % estimate gravity field (TBD)
        hc = Hc(:, 2:end);
        C = null(Hp');
       
        total_inf = total_inf + ((C'*hc)' * ((C'*R0*C)\(C'*hc)));
        
         if mod(k,Nbatch) == 0
            flag = 1;     % trigger
            flag_count = flag_count + 1;
         else
            flag = 0;     % reset
         end
        if(flag)
            n_solve = fix((rank(total_inf) + 4)^(1/2) - 1);

            Nt = k - kmin + 1;

            % calibrated signal
            Y_cal = Y(:, kmin:k) - X(7:end, kmin:k);

            [X_CS, P0c_new, NSM_err] = estimate_grav_field(n_solve, time(kmin:k), ...
                Nt, Y_cal, X(1:3, kmin:k), X(4:6, kmin:k), ...
                planetParams, poleParams, Cnm, Snm, P0c_grav, R0,...
                Px(kmin:k, :));
            delta_grav = (0.3).*abs(P0c_grav(2:end, 2:end) - P0c_new);
            
            % update uncertainty grav. field + coefficient value
            [Cnm, Snm, P0c, P0c_grav] = update_coeff(n_solve, planetParams(3), ...
                X_CS, Cnm, Snm, delta_grav + P0c_new, P0c, P0c_grav);
            
            % update uncertainty navigation
            p = P(1:Ns,1:Ns);
            P = [p, Pxc;Pcx, P0c];

            % display NSM filter error
            fprintf('   Gravity estimation number %d. Initial error: %.2e. Final error: %.2e.\n', ...
            flag_count, NSM_err(1), NSM_err(end));

            kmin = k+1;
        end
        % compute coefficient error
        Xc     = mat2list(Cnm, Snm, Ncoeff, Nsoeff);
        errC   = abs(Xc - Xc_t);
        sigmaC = diag(sqrt(P0c_grav));

        % Update plots
        [errNorm_pos, errNorm_vel, ...
        errNorm_bias, std_pos, std_vel, std_bias] = ...
        update_plots(time, k, state_true, X0, P, h, ax, errNorm_pos,...
        errNorm_vel, errNorm_bias, std_pos, std_vel, std_bias, ...
        errC, sigmaC, planetParams(3));
    
        % stop if figure is closed
        if ~ishandle(fig)
            disp('Figure closed, stopping simulation.');
            break;
        end
    end
end


%% Auxiliary functions
function [Cnm, Snm, P0c, P0c_grav] = update_coeff(n_solve, n_max, X_CS, Cnm, Snm,...
    P0, P0c, P0c_grav)

    [Nc, Ns, Ncs]  = count_num_coeff(n_solve); 
    [Nc_org, ~, ~]  = count_num_coeff(n_max); 
    [Cnm_new, Snm_new] = list2mat(n_solve, Nc, Ns, X_CS);

    Cnm(1:n_solve+1,1:n_solve+1) = Cnm_new;
    Snm(1:n_solve+1,1:n_solve+1) = Snm_new;
    
    P0cc = diag(P0(1:Nc-1, 1:Nc-1)); P0cs = P0(1:Nc-1, Nc:Ncs-1);
    P0sc = P0(Nc:Ncs-1,1:Nc-1); P0ss = diag(P0(Nc:Ncs-1, Nc:Ncs-1));

    P0c(2:Nc, 2:Nc) = diag(P0cc);
    P0c(2:Nc, Nc_org+1:Ns+Nc_org) = P0cs.*0;
    P0c(Nc_org+1:Ns+Nc_org, 2:Nc) = P0sc.*0;
    P0c(Nc_org+1:Ns+Nc_org, Nc_org+1:Ns+Nc_org) = diag(P0ss);

    P0c_grav(2:Nc, 2:Nc) = P0(1:Nc-1, 1:Nc-1);
    P0c_grav(2:Nc, Nc_org+1:Ns+Nc_org) = P0cs;
    P0c_grav(Nc_org+1:Ns+Nc_org, 2:Nc) = P0sc;
    P0c_grav(Nc_org+1:Ns+Nc_org, Nc_org+1:Ns+Nc_org) = P0(Nc:Ncs-1, Nc:Ncs-1);

    P0c_grav = diag(diag(P0c_grav));
end