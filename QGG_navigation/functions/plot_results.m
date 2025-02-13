function plot_results(t, state_true, X, P, Pc, Xhat, pref, posf, planetParams, ...
    count, considerCov, augmented_st, posE, posM, system)
    %%                    PLOT RESULTS FUNCTION
    % Description: Plot the errors + uncertainty in the estimation.
    % Author: Sergio Coll Ibars
    % Date: 03/29/2024
    
    % tselect and format time
    if(system == "EPHEM")
        jd = 2451545 + t / planetParams(3) / 86400;
        humanReadableTime = datetime(jd, 'ConvertFrom', ...
            'juliandate');
        humanReadableTime.Format = 'MMM dd, yyyy';
        date_init = string(humanReadableTime(1));
        date_end  = string(humanReadableTime(end));
        humanReadableTime.Format = 'MMM dd';

        time = humanReadableTime';
        xlb = "date";
        tt  = ' from ' + date_init + ' - ' + date_end;
    else
        time = t'/planetParams(3)/86400;
        xlb = "days";
        tt = '.';
    end
    At = (t(2) - t(1)); % [-]

    Ns = length(X(:, 1));
    
    % compute Earth and Moon motion
    rE = posE;
    rM = posM;
    vM = [diff(posM(1, :));diff(posM(2, :)); diff(posM(3, :))]/At;
    vM = [zeros(3, 1), vM];

    % plot options
    lw = 2;
    color1 = [204, 0, 204]./256;     % violet
    color2 = "#FF0000";              % red
    color3 = "k";                    % black
    set(0,'defaultAxesFontSize',16);

    % plot real trajectory and reconstructed
    figure()
    plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
        'LineWidth', lw, 'Color', color3)
    hold on;
    plot3(X(1, :), X(2, :), X(3, :), 'LineWidth', lw, 'Color', color2)
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title("true vs reconstructed trajectory")
    legend("true", "reconstructed")
    axis equal;
    grid on;
    
     % plot real trajectory and reconstructed
    figure()
    plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
        'LineWidth', lw, 'Color', color3)
    hold all;
    plot3(rE(1, :), rE(2, :), rE(3, :), 'LineWidth', lw, 'Color', color2, ...
        'LineStyle','--')
    plot3(rM(1, :), rM(2, :), rM(3, :), 'LineWidth', lw, 'Color', color1, ...
        'LineStyle','--')
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title("3BP trajectories. Inertial frame")
    legend("S/C", "Earth", "Moon")
    axis equal;
    grid on;

    % convert to RNT frame
    r = X(1:3, :) - rM;
    v = X(4:6, :) - vM;
    [~, ~, P_RTN] = convert2RNT(r, v, t, P, Ns);

    % plot error
    err = state_true(:, 1:Ns)' - X(1:Ns, :);
    sigma = sqrt(diag(Pc));
<<<<<<< HEAD
    eigP = zeros(Ns, length(t));
=======
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
    cov = zeros(Ns, length(t)); cov2 = cov; cov2_RTN = cov;
    cov2_tot = zeros(6, length(t));
    cov2_E = zeros(6, length(t));
    cov2_M = zeros(6, length(t));
    Np = sqrt(length(P(1, :)));
    for j = 1:length(t)
        p = reshape(P(j, :), [Np,Np]); p_RTN = reshape(P_RTN(j, :), [Np,Np]);
        a = sqrt(diag(p)); a_RTN = sqrt(diag(p_RTN));
<<<<<<< HEAD
        eigP(: , j) = eig(p);
=======
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
        cov(:, j) = a(1:Ns);
        cov2(:, j) = a(1:Ns).^2;
        cov2_RTN(:, j) = a_RTN(1:Ns).^2;
        if(considerCov)
            pxc_tot = p(1:6, 7:end);
            Sxc_tot = pxc_tot/(Pc);
            n = length(pxc_tot(1, :));
            Sxc_E = Sxc_tot(:, 1:n/2);
            Sxc_M = Sxc_tot(:, n/2+1:end);
            cov2_tot(:, j) = sum(Sxc_tot * diag(sigma), 2);
            cov2_E(:, j) = sum(Sxc_E * diag(sigma(1:n/2)), 2);
            cov2_M(:, j) = sum(Sxc_M * diag(sigma(n/2+1:end)), 2);
        end
    end

    figure()
    subplot(1, 2, 1)
    d = sqrt(sum(cov2(1:3, :), 1));
    scale = planetParams(2)./1000;
    semilogy(time, vecnorm(err(1:3, :)).* scale, 'LineWidth', lw, 'Color', color1)
    hold on;
    semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
    grid on;
    ylabel('[Km]')
    title('Position norm')

    subplot(1, 2, 2)
    d = sqrt(sum(cov2(4:6, :), 1));
    scale = planetParams(3) * planetParams(2);
    semilogy(time, vecnorm(err(4:6, :)).* scale, 'LineWidth', lw, 'Color', color1)
    hold on;
    semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
    grid on;
    ylabel('[m/s]')
    title('Velocity norm')
    legend('error', '3 \sigma')
    sgtitle('State error vector norm and 3 \sigma bound')

<<<<<<< HEAD
    figure()
    scale1 = planetParams(2);                   % [m]
    scale2 = planetParams(2)*planetParams(3);   % [m/s]
    scale = [scale1, scale1, scale1, scale2, scale2, scale2];
    index = [1, 3, 5, 2, 4, 6];
    for j = 1:6
        subplot(3, 2, index(j))
        semilogy(time, sqrt(eigP(j, :))*scale(j), 'LineWidth', 2, 'Color', 'b')
        grid on;
    end
    xlabel('Date')
    sgtitle('Covariance eigen-value over time')

=======
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074

    figure()
    index = [1, 3, 5];
    scale = planetParams(2)./1000;
    for j = 1:3
        subplot(3, 2, index(j))
        upper_bound = +3*cov(j, :) * scale;
        lower_bound = -3*cov(j, :) * scale;
        plot(time, err(j, :) * scale, 'LineWidth', lw, 'Color', color1)
        hold on;
        fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
            color1, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        xlabel(xlb)
        ylabel('error R_' + string(j) + '[km]')
        
        str = {'RMS = ' + string(rms(err(j, :) * scale))};
        legend(str,'Location', 'best');
    end

    index = [2, 4, 6];
    scale = planetParams(3) * planetParams(2);
    for j = 1:3
        subplot(3, 2, index(j))
        upper_bound = +3*cov(j+3, :) * scale;
        lower_bound = -3*cov(j+3, :) * scale;
        plot(time, err(j+3, :) * scale, 'LineWidth', lw, 'Color', color1)
        hold on;
        fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
            color1, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        xlabel(xlb)
        ylabel('error V_' + string(j) + '[m/s]')

        str = {'RMS = ' + string(rms(err(j+3, :) * scale))};
        legend(str,'Location', 'best');
    end
    sgtitle("State error + 3\sigma bounds " + tt)
    
    figure()
    index = [1, 3, 5];
    scale = planetParams(2)./1000;
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, abs(err(j, :)) * scale, 'LineWidth', lw, 'Color', color1)
        xlabel(xlb)
        ylabel('error R_' + string(j) + '[km]')
       
        
        str = {'RMS = ' + string(rms(err(j, :) * scale))};
        legend(str,'Location', 'best');
    end

    index = [2, 4, 6];
    scale = planetParams(3) * planetParams(2);
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, abs(err(j+3, :)) * scale, 'LineWidth', lw, 'Color', color1)
        xlabel(xlb)
        ylabel('error V_' + string(j) + '[m/s]')

        str = {'RMS = ' + string(rms(err(j+3, :) * scale))};
        legend(str,'Location', 'best');
    end
    sgtitle("State error" + tt)

    % plot state correction
    figure()
    index = [1, 3, 5];
    scale = planetParams(2)./1000;
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, Xhat(j, :) * scale, 'LineWidth', lw, 'Color', color2)
        xlabel(xlb)
        ylabel('corr R_' + string(j) + '[Km]')
    end

    index = [2, 4, 6];
    scale = planetParams(3) * planetParams(2);
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, Xhat(j+3, :)*scale, 'LineWidth', lw, 'Color', color2)
        ylabel('corr V_' + string(j) + '[m/s]')
    end
    sgtitle('correction Xhat vector')


    % plot  trace uncertainty. Inertial frame
    figure()
    subplot(1, 2, 1)
    semilogy(time, sum(cov(1:3, :)) * planetParams(2), ...
        'LineWidth', lw, 'Color', 'g')
    xlabel(xlb)
    ylabel('[m]')
    title('position')

    subplot(1, 2, 2)
    semilogy(time, sum(cov(4:6, :)) * planetParams(2) * planetParams(3),...
        'LineWidth', lw, 'Color', 'g')
    xlabel(xlb)
    ylabel('[m/s]')
    title('velocity')
    sgtitle('Covariance trace over time. Inertial frame.  ' + tt)
    grid on;

    % plot  trace uncertainty. RTN frame
    figure()
    index = [1, 3, 5, 2, 4, 6];
<<<<<<< HEAD
    ttitle = ["R", "T", "N", "V_R", "V_T", "V_N"];
=======
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
    for k = 1:6
        subplot(3, 2, index(k));
        semilogy(time, sqrt(cov2_RTN(k, :)), ...
            'LineWidth', lw, 'Color', 'g')
        xlabel(xlb)
        ylabel('[-]')
<<<<<<< HEAD
        title(ttitle(k));
    end
    sgtitle('Covariance evolution over time. RTN Moon centered frame')
=======
    end
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074


    % plot 1 sigma state error due to consider parameters
    figure()
    index = [1, 3, 5, 2, 4, 6];
    scale = planetParams(2)./1000;
    ttitle = {'\Delta_x [km]', '\Delta_y [km]', '\Delta_z [km]', ...
        '\Delta_{vx} [m/s]', '\Delta_{vy} [m/s]',  '\Delta_{vz} [m/s]'};
    for j =1:3
        subplot(3, 2, index(j))
        plot(time, cov2_tot(j, :)*scale, time, cov2_E(j, :)*scale, time, ...
            cov2_M(j, :)*scale, 'LineWidth', lw)
        title(ttitle(j))
        xlabel(xlb)
        if j ==1, legend('Total', 'Earth coeff', 'Moon coeff'), end
    end
    scale = planetParams(3) * planetParams(2);
    for j =4:6
        subplot(3, 2, index(j))
        plot(time, cov2_tot(j, :)*scale, time, cov2_E(j, :)*scale, time, ...
            cov2_M(j, :)*scale, 'LineWidth', lw)
        title(ttitle(j))
        xlabel(xlb)
    end
    sgtitle('Errors in state estimate due to 1\sigma error in consider parameters');

    % plot both trajectories
    figure()
    index = [1, 3, 5];
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, state_true(:, j), 'LineWidth', lw, 'Color', color1)
        hold on;
        plot(t, X(j, :), 'LineWidth', lw, 'Color', color2)
        xlabel(xlb)
        ylabel('R_' + string(j) + '[-]')
        if(j ==1)
            legend('true', 'reconstructed')
        end
    end

    index = [2, 4, 6];
    for j = 1:3
        subplot(3, 2, index(j))
        plot(time, state_true(:, j+3), 'LineWidth', lw, 'Color', color1)
        hold on;
        plot(t, X(j+3, :), 'LineWidth', lw, 'Color', color2)
        xlabel(xlb)
        ylabel('V_' + string(j) + ' [-]')
    end
    sgtitle("True and reconstructed trajectory " + tt)

    % plot prefit & posfit per iter
    ncols = floor(sqrt(count));
    nrows = ceil(count/ncols);
    measDim = planetParams(3)^2;      % [Etvos]
    if(augmented_st), Ns  = 6; else, Ns = 7; end
    figure()
    for j = 1:count
        maxInd = Ns * j;
        minInd = maxInd - (Ns - 1);

        subplot(nrows, ncols, j)
        plot(time, pref(minInd:maxInd, :) * measDim, '.')
        title('Iter ' + string(j))
        ylabel('[1/s^2]')
    end
    sgtitle('Prefit per iteration')
    
    figure()
    for j = 1:count
        maxInd = Ns * j;
        minInd = maxInd - (Ns-1);

        subplot(nrows, ncols, j)
        plot(time, posf(minInd:maxInd, :) * measDim, '.')
        title('Iter ' + string(j))
    end
    sgtitle('Posfit per iteration')

    maxInd = count * Ns;
    minInd = maxInd - (Ns-1);
    figure()
    subplot(1, 2, 1)
    plot(time, pref(minInd:maxInd, :) * measDim, '.')
    xlabel(xlb)
    title('prefit iter = '  + string(count))
    subplot(1, 2, 2)
    plot(time, posf(minInd:maxInd, :) * measDim, '.')
    xlabel(xlb)
    title('postfit iter = '  + string(count))

% %     % plot unmodeled acceleration
% %     scale = planetParams(2) * planetParams(3);  % [m/s]
% %     if(length(X(:, j)) > 6)
% %         figure()
% %         plot(time, X(7:end, :).*scale, 'LineWidth', 2)
% %         title('Unmodeled acceleration over time')
% %     end

    % plot augmented state estimation
    if(augmented_st)
        figure()
        semilogy(time, X(7:end, :), 'LineWidth', 2)
        hold on;
        semilogy(time, abs(err(7, :)), 'LineWidth', 2, 'Color', 'r')
        semilogy(time, 3*cov(7, :), time, -3*cov(7, :), 'LineWidth', 2, 'Color',...
            'k', 'LineStyle', '--')
        title('SRP estimation, \eta factor')
        xlabel(xlb)
        ylabel('[-]')
    end
end

