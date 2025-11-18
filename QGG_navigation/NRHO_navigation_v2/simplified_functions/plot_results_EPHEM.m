function [] = plot_results_EPHEM(t, state_true, X, P, pref, posf,...
    planetParams, count, r_sc_moon)
    % Plot the estimation results in the EPHEMERIDES case.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state number 
    Ns = 12;    % number of states
    Nm = 6;     % number of measurements

    % get apo-apsis and peri-apsis passages
    idx_min = islocalmin(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations
    idx_max = islocalmax(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations

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
    hold all;
    semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
    if(sum(idx_min)~=0)
        xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
            'LineStyle', '--')
    end
    if(sum(idx_max)~=0)
        xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
            'LineStyle', '--')
    end
    grid on;
    ylabel('[Km]')
    title('Position norm')

    subplot(1, 2, 2)
    d = sqrt(sum(cov2(4:6, :), 1));
    scale = planetParams(3) * planetParams(2);
    semilogy(time, vecnorm(err(4:6, :)).* scale, 'LineWidth', lw, 'Color', color1)
    hold all;
    semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k');
    if(sum(idx_min)~=0)
        xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
            'LineStyle', '--')
    end
    if(sum(idx_max)~=0)
        xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
            'LineStyle', '--')
    end
    grid on;
    ylabel('[m/s]')
    title('Velocity norm')
    legend('error', '3 \sigma')
    sgtitle('State error vector norm and 3 \sigma bound')

    %% plot gradiometer bias + 3 sigma
    figure()
    ttlabel = ["xx", "xy", "xz", "yy", "yz", "zz"];
    for k = 1:6
        subplot(3, 2, k)
        d  = cov(6+k, :);
        scale = 1E3;
        plot(time, err(6+k, :).* scale, 'LineWidth', 1.2, ...
            'Color', color2)
        hold all;
        plot(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
        plot(time, -3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
        grid on;
        ylabel('[Eotvos]')
        title('B_{' + ttlabel(k) + '}');
    end
    sgtitle('Gradiometer bias estimation error');


    figure()
    ttlabel = ["xx", "xy", "xz", "yy", "yz", "zz"];
    for k = 1:6
        subplot(3, 2, k)
        scale = 1E3;
        plot(time, X(6+k, :).* scale, 'LineWidth', 1.2, 'Color', 'r')
        hold all;
        plot(time,state_true(:, 6+k).* scale, 'LineWidth', lw, ...
            'Color', 'k')
        grid on;
        ylabel('[Eotvos]')
        title('B_{' + ttlabel(k) + '}');
    end
    sgtitle('True + predicted bias error');

    %% prefit + postfit per iteration
    ncols = floor(sqrt(count));
    nrows = ceil(count/ncols);
    measDim = 1;      % [Eotvos]
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

    %% plot condition number for the uncertainty over time
    condNum = ones(1, length(t)) * NaN;
    eigenVals = ones(12, length(t)) * NaN;
    for k =1:length(t)
        p = reshape(P(k, :), [Np,Np]); 
        condNum(k) = cond(p);
        eigenVals(:, k) = eig(p);
    end
    figure()
    semilogy(time, condNum, 'LineWidth', 2)
    grid on;
    xlabel('Time');
    title('Condition number for formal uncertainty');

    figure()
    for k = 1:12
        subplot(3, 4, k)
        semilogy(time, eigenVals(k, :), 'LineWidth', 2)
        grid on;
    end
    sgtitle('Eigen-values covariance matrix')
end

