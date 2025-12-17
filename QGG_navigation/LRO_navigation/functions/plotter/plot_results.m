function [] = plot_results(time, state_true, bias_true, X_EKF, P_EKF, mask)
    %% PLOT FILTER RESULTS
    
    % get time and states numbers
    Nt  = length(time);
    Nx  = length(X_EKF(:, 1));

    % convert time to UTC
    utc = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

    % get estimation formal error
    sigma = nan(Nx, Nt);
    for k = 1:Nt
        p = reshape(P_EKF(k, :), [Nx, Nx]);
        sigma(:, k) = sqrt(diag(p));
    end

    % plot error position [m]
    figure();
    errorP = state_true(:, 1:3)' - X_EKF(1:3, :);
    tt = ["X", "Y", "Z"];
    for k = 1:3
        subplot(1, 3, k);
        plot(tUTC, errorP(k, :), 'LineWidth', 2); hold all;
        plot(tUTC, +3.*sigma(k, :), 'LineWidth', 2, 'Color', 'k');
        plot(tUTC, -3.*sigma(k, :), 'LineWidth', 2, 'Color', 'k');
        grid on;
        ylabel('[m]'); title(tt(k));
    end
    sgtitle('position error')

    % plot error velocity [m/s]
    figure();
    errorV = state_true(:, 4:6)' - X_EKF(4:6, :);
    tt = ["V_x", "V_y", "V_z"];
    for k = 1:3
        subplot(1, 3, k);
        plot(tUTC, errorV(k, :), 'LineWidth', 2); hold all;
        plot(tUTC, +3.*sigma(k+3, :), 'LineWidth', 2, 'Color', 'k');
        plot(tUTC, -3.*sigma(k+3, :), 'LineWidth', 2, 'Color', 'k');
        grid on;
        ylabel('[m/s]'); title(tt(k));
    end
    sgtitle('velocity error')

    % plot bias error
    figure();
    tt   = ["xx", "xy", "xz", "yy", "yz", "zz"];
    errB = bias_true - X_EKF(7:end, :); 
    for k = 1:6
        subplot(2, 3, k);
        if(mask(k) == 1)
            plot(tUTC, errB(k, :), 'LineWidth', 2); hold all;
            plot(tUTC, +3.*sigma(6+k, :), 'LineWidth', 2, 'Color','k');
            plot(tUTC, -3.*sigma(6+k, :), 'LineWidth', 2, 'Color','k');
        end
        title(tt(k)); ylabel('[mE]'); grid on;
    end
    sgtitle('Bias estimation error')
end

