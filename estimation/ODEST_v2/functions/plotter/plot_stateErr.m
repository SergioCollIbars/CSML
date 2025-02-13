function [] = plot_stateErr(stateErr, stateSigma, time, count)
    %PLOT_STATEERR plot position and velocity errors at each 
    % iteration number.
    
    % set time to days
    time = time./86400; % [days]

    % sobplot rows and cols number
    nrows = floor(sqrt(count+1));
    ncols = ceil((count+1)/nrows);

    % lines colors
    c = ['m', 'b', 'g', 'm', 'b', 'g'];

    figure()
    val = ones(6, length(stateErr(1, :, 1)));
    sigma = ones(6, length(stateErr(1, :, 1)));
    for j = 1:count + 1
        val(:, :) = stateErr(:, :, j);
        sigma(:, :) = stateSigma(:, :, j);
        
        subplot(nrows, ncols, j);
        for k = 1:3
            plot(time, val(k, :), 'LineWidth', 2, 'Color', c(k))
            hold on;
            if(j > 1)
            plot(time, 3.*sigma(k, :), time, -3.*sigma(k, :), 'LineWidth', 2, ...
                'LineStyle', '--', 'Color', c(k))
            end
        end
                    xlabel('Time [days]')
            ylabel('[m]')
        title('Iteration = ' + string (j - 1))
        if(j == 1), legend('x', 'y', 'z'), end
    end
    sgtitle('Position errors')

    figure()
    val = ones(6, length(stateErr(1, :, 1)));
    sigma = ones(6, length(stateErr(1, :, 1)));
    for j = 1:count + 1
        val(:, :) = stateErr(:, :, j);
        sigma(:, :) = stateSigma(:, :, j);

        subplot(nrows, ncols, j);
        for k = 4:6
            plot(time, val(k, :), 'LineWidth', 2, 'Color', c(k))
            hold on;
            if(j > 1)
            plot(time, 3.*sigma(k, :), time, -3.*sigma(k, :), 'LineWidth', 2, ...
                'LineStyle', '--', 'Color', c(k))
            end

        end
                    xlabel('Time [days]')
            ylabel('[m/s]')
        title('Iteration = ' + string (j - 1))
        if(j == 1), legend('v_x', 'v_y', 'v_z'), end
    end
    sgtitle('Velocity errors')
end

