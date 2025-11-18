function [] = plot_results_EPHEM(t, state_true, X, P, pref, posf,...
    planetParams, count)
    % Plot the estimation results in the EPHEMERIDES case.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state number 
    Ns = 16;    % number of states
    Nm = 6;     % number of measurements

    % plot options
    lw = 2;
    color1 = [204, 0, 204]./256;     % violet
    color2 = "#FF0000";              % red
    color3 = "k";                    % black
    set(0,'defaultAxesFontSize',16);

    % convert to date time
    jd = 2451545 + t / planetParams(3) / 86400;
    humanReadableTime = datetime(jd, 'ConvertFrom', ...
        'juliandate');
    humanReadableTime.Format = 'MMM dd, yyyy';
    date_init = string(humanReadableTime(1));
    date_end  = string(humanReadableTime(end));
    humanReadableTime.Format = 'MMM dd';

    time = humanReadableTime';
    
    %% parse error and covariance values
    err = state_true(:, 1:Ns)' - X(1:Ns, :);
    cov = zeros(Ns, length(t)); cov2 = cov;
    Np = sqrt(length(P(1, :)));
    for j = 1:length(t)
        p = reshape(P(j, :), [Np,Np]); 
        a = sqrt(diag(p));

        cov(:, j)  = a(1:Ns);
        cov2(:, j) = a(1:Ns).^2;
    end 

    %% plot real trajectory and reconstructed
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

    %% Plot state error + 3 sigma
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

    %% plot SRP scale factor + 3 sigma
    figure()
    semilogy(time, X(7, :), 'LineWidth', 2)
    hold on;
    semilogy(time, abs(err(7, :)), 'LineWidth', 2, 'Color', color2)
    semilogy(time, 3*cov(7, :), time, -3*cov(7, :), 'LineWidth', 2, 'Color',...
            'k', 'LineStyle', '--')
    title('SRP estimation, \eta factor')
    ylabel('[-]')


    %% plot gradiometer bias + 3 sigma
    figure()
    ttlabel = ["xx", "xy", "xz", "yy", "yz", "zz"];
    for k = 1:6
        subplot(3, 2, k)
        d  = cov(7+k, :);
        scale = planetParams(3)^2/1E-9;
        semilogy(time, abs(err(7+k, :)).* scale, 'LineWidth', lw, ...
            'Color', color2)
        hold on;
        semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
        grid on;
        ylabel('[Eotvos]')
        title('B_{' + ttlabel(k) + '}');
    end
    sgtitle('Gradiometer bias estimation error');

    %% plot accelerometer bias + 3 sigma
    figure()
    ttlabel = ["x", "y", "z"];
    for k = 1:3
        subplot(1, 3, k)
        d = cov(13+k, :);
        scale = planetParams(2) * planetParams(3);
        semilogy(time, abs(err(13+k, :)).* scale, 'LineWidth', lw, ...
            'Color', color2)
        hold on;
        semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
        grid on;
        ylabel('[m/s^2]')
        title('B_{' + ttlabel(k) + '}');
    end
    sgtitle('Accelerometer bias estimation error');

    %% prefit + postfit per iteration
    ncols = floor(sqrt(count));
    nrows = ceil(count/ncols);
    measDim = planetParams(3)^2/1E-9;      % [Eotvos]
    figure()
    for j = 1:count
        maxInd = Nm * j;
        minInd = maxInd - (Nm - 1);

        subplot(nrows, ncols, j)
        plot(time, pref(minInd:maxInd, :) * measDim, '.')
        title('Iter ' + string(j))
        ylabel('[1/s^2]')
    end
    sgtitle('Prefit per iteration');
    
    figure()
    for j = 1:count
        maxInd = Nm * j;
        minInd = maxInd - (Nm-1);

        subplot(nrows, ncols, j)
        plot(time, posf(minInd:maxInd, :) * measDim, '.')
        title('Iter ' + string(j))
    end
    sgtitle('Posfit per iteration');
end

