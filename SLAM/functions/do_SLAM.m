function [X, Px, Xc, P0c_grav] = do_SLAM(time, state_true, planetParams, ...
    poleParams, instrumentParams, BN_mat, Cnm_t, Snm_t, Y)
    
    % load filter parameters
    [R0, P0,  P0c_grav, Q0, Qb, delta_state0, Cnm, Snm, ~] = ...
        loadFilterParams(planetParams, instrumentParams, Cnm_t, Snm_t);
    
    % initialize filter
% %     [P0, state0] = initialize_filter(time, ...
% %         delta_state0+state_true(1, 1:12)', Y, R0,...
% %         P0, 20, planetParams, poleParams, Cnm, Snm, BN_mat);
% % 
% %     err0 = state0 - state_true(1, 1:12)';
% %     s0   = 3 * sqrt(diag(P0));

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

    % Create figure and graphics objects
    fig = figure('Name','EKF State Error & Covariance',...
        'NumberTitle','off');
    subplot(2,1,1); hold on; grid on;
    hErr_pos = plot(nan, nan, '-');  % handle for error curve
    hCov_pos = plot(nan, nan, '-');  % handle for error curve
    xlabel('Time [hr]');
    ylabel('[m]');
    title('Position error Norm');
    xlim([0, time(end)./3600]);
    ax1 = gca;
    
    subplot(2,1,2); hold on; grid on;
    hErr_vel = semilogy(nan, nan, '-');  % handle for error curve
    hCov_vel = semilogy(nan, nan, '-');  % handle for error curve
    xlabel('Time [hr]');
    ylabel('[mm/s]');
    title('Velocity error Norm');
    xlim([0, time(end)./3600]);
    ax2 = gca;

    figure('Name','Bias errors','NumberTitle','off');
    h = gobjects(1, 6);
    c = gobjects(1, 6);
    ax= gobjects(1, 6);
    
    % --- Row 1 ---
    for f = 1:6
        ax(f) = subplot(3,2,f); hold on; grid on;
        h(f) = plot(nan, nan, 'b-', 'LineWidth', 1.5);
        c(f) = plot(nan, nan, 'b-', 'LineWidth', 1.5);
        title('Bias ' + string(f)); xlabel('Time [hr]'); ylabel('[mE]');
        xlim([0, time(end)./3600]);
    end

    figure('Name','EKF Grav estimation error',...
        'NumberTitle','off');
    hold on; grid on;
    semilogy(1:Ncs-1, abs(Xc_t(2:end)), 'LineWidth', 2, 'Color','b');
    hErr_c = semilogy(nan, nan, '-');  % handle for error curve
    hCov_c = semilogy(nan, nan, '-');  % handle for error curve
    xlabel('Coeff');
    ylabel('[-]');
    title('Gravity coeff error');
    ax3 = gca;
    set(ax3, 'YScale','log');
    
    errC   = abs(Xc_0 - Xc_t);
    sigmaC = sqrt(diag(P0c_grav));

    set(hErr_c, 'XData', 1:Ncs-1, 'YData',...
        errC(2:end), 'LineWidth', 1.5, 'Color', 'r');
    set(hCov_c, 'XData', 1:Ncs-1, 'YData', 3.*sigmaC(2:end), ...
        'LineWidth', 2, 'Color', 'k');
    
    total_inf = inv(P0c_grav(2:end, 2:end)); kmin = 1;
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
        
         if mod(k,1000) == 0
            flag = 1;     % trigger
         else
            flag = 0;     % reset
         end
        if(flag)
            n_solve = fix((rank(total_inf) + 4)^(1/2) - 1);

            Nt = k - kmin + 1;
            sigmaPos  = max((std_pos(:,kmin:k)./3)');
            sigmaBias = max((std_bias(:,kmin:k)./3)');

            % calibrated signal
            Y_cal = Y(:, kmin:k) - X(7:end, kmin:k);

            [X_CS, P0c_new] = estimate_grav_field(n_solve, time(kmin:k), ...
                Nt, Y_cal, X(1:3, kmin:k), X(4:6, kmin:k), ...
                planetParams, poleParams, Cnm, Snm, P0c_grav, R0.*(100),...
                sigmaPos, sigmaBias);

            [Cnm, Snm, P0c, P0c_grav] = update_coeff(n_solve, planetParams(3), ...
                X_CS, Cnm, Snm, P0c_new, P0c, P0c_grav);

            % % p = diag(P);
            % % p(7:end) = p(7:end).*(1.3);
            % % P = [diag(p(1:Ns)), Pxc;Pcx, P0c.*1];
            p = P(1:Ns,1:Ns);
            P = [p, Pxc;Pcx, P0c];

            % compute coefficient error
            Xc     = mat2list(Cnm, Snm, Ncoeff, Nsoeff);
            errC   = abs(Xc - Xc_t);
            sigmaC = diag(sqrt(P0c_grav));

            set(hErr_c, 'XData', 1:Ncs-1, 'YData',...
            errC(2:end), 'LineWidth', 1.5, 'Color', 'r');
            set(hCov_c, 'XData', 1:Ncs-1, 'YData', 3.*sigmaC(2:end), ...
            'LineWidth', 2, 'Color', 'k');

            kmin = k+1;
        end

        % Update plots
        [errNorm_pos, ...
            errNorm_vel, errNorm_bias, std_pos, std_vel, ...
            std_bias]= ...
            update_openPlots(time, k, state_true, X0, P(1:Ns, 1:Ns), ...
            hCov_pos, hCov_vel, hErr_pos, hErr_vel, h, c, errNorm_pos, ...
            errNorm_vel, errNorm_bias, std_pos, std_vel, std_bias,...
            ax1, ax2, ax);
    
        % stop if figure is closed
        if ~ishandle(fig)
            disp('Figure closed, stopping simulation.');
            break;
        end
    end

% %     % save gravity field uncertainty
% %     disp('Save Gravity data ...');
% %     save("data/gravField_cov.mat", 'P0c_grav');
% %     save("data/gravField_coeff.mat", 'Xc');
end


%% Auxiliry functions
function [errNorm_pos, errNorm_vel, errNorm_bias, ...
    std_pos, std_vel, std_bias] = update_openPlots(time, k, state_true, X0, P,...
    hCov_pos, hCov_vel, hErr_pos, hErr_vel, h, c, errNorm_pos, errNorm_vel, ...
    errNorm_bias, std_pos, std_vel, std_bias, ax1, ax2, ax)
        e_k            = state_true(k, 1:6)'  - X0(1:6);
        e_b            = state_true(k, 7:12)' - X0(7:12); 
        errNorm_pos(k) = norm(e_k(1:3)); errNorm_vel(k) = norm(e_k(4:6));
        errNorm_bias(:, k) = abs(e_b);
        std = diag(sqrt(P));
        std_pos(k)     = 3*norm(std(1:3)); std_vel(k)   = 3*norm(std(4:6));
        std_bias(:, k) = 3.*std(7:12);

        % update X and Y data of the existing line objects
        set(hErr_pos, 'XData', time(1:k)./3600, 'YData',...
            errNorm_pos(1:k), 'LineWidth', 1.5, 'Color', 'r');
        set(hCov_pos, 'XData', time(1:k)./3600, 'YData', std_pos(1:k), ...
            'LineWidth', 2, 'Color', 'k');
        set(ax1, 'YLim', [0, 5*std_pos(k)]);

        set(hErr_vel, 'XData', time(1:k)./3600, 'YData',...
            errNorm_vel(1:k).*1E3, 'LineWidth', 1.5, 'Color', 'r');
        set(hCov_vel, 'XData', time(1:k)./3600, 'YData', std_vel(1:k).*1E3, ...
            'LineWidth', 2, 'Color', 'k');
        set(ax2, 'YLim', [0, 5*std_vel(k)*1E3]);

        for f = 1:6
            set(h(f), 'XData', time(1:k)./3600, 'YData',...
            errNorm_bias(f, 1:k), 'LineWidth', 1.5, 'Color', 'r');
            set(c(f), 'XData', time(1:k)./3600, 'YData', std_bias(f, 1:k), ...
            'LineWidth', 2, 'Color', 'k');
            set(ax(f), 'YLim', [0, 5*std_bias(f, k)]);
        end
    
        drawnow limitrate;  % update the figure without killing performance
end

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
end