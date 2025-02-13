function [] = plot_coeffErr(n_max,RMSerr, RMSsigma, coeffErr, count)
    % PLOT_COEFFERR plot SH coefficient errors. RMS and SRN
    
    % disable warnings
    warning('off')

    % sobplot rows and cols number
    nrows = floor(sqrt(count+1));
    ncols = ceil((count+1)/nrows);
    [Nc, Ns, ~] = count_num_coeff(n_max);

    figure()
    for j = 1:count + 1
        val   = RMSerr(j, 2:end);
        sigma = RMSsigma(j, 2:end);

        subplot(nrows, ncols, j);
        errorbar(2:n_max, val, sigma, 'LineWidth', 2, 'Marker','o', ...
            'LineStyle','--', 'Color', 'b')
        set(gca,'YScale','log')
        xlabel('SH degree [-]')
        ylabel('RMS error [-]')
        title('Iteration = ' + string (j - 1))
    end
    sgtitle('RMS degree error SH coefficients')

    % enable warnings
    warning('on')

    % compute coefficient errors
    figure()
    plot(1:Nc, coeffErr(count+1, 1:Nc).*100, 'Marker','o', 'Color', 'b', ...
        'MarkerFaceColor', 'b')
    grid on;
    xlabel('Degree')
    ylabel('[%]')
    title('Percentage error for C_{nm} SH coefficients')
    ylim([0, 100])
    
    figure()
    plot(1:Ns, coeffErr(count+1, Nc+1:Nc+Ns).*100, 'Marker','o', 'Color', 'b', ...
        'MarkerFaceColor', 'b')
    grid on;
    xlabel('Degree')
    ylabel('[%]')
    title('Percentage error for S_{nm} SH coefficients')
    ylim([0, 100])
end

