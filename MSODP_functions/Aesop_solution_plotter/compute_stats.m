function [] = compute_stats(rel_err, tt)
    % stats
    valid = rel_err(~isnan(rel_err));
    
    % Total number of valid entries
    N = numel(valid);
    
    % Count > 1
    count_gt1 = sum(valid > 1);
    
    % Count < 1
    count_lt1 = sum(valid < 1);
    
    % Percentages
    perc_gt1 = 100 * count_gt1 / N;
    perc_lt1 = 100 * count_lt1 / N;
    
    % Display
    fprintf('   Percentage > 1: %.2f%%\n', perc_gt1);
    fprintf('   Percentage < 1: %.2f%%\n', perc_lt1);
    
    % plot Histogram
    figure();
    histogram(log10(valid), 100); hold all
    MEDIAN = median(log10(valid));
    ylabel('# of Events');
    xlabel('Relative error');
    xticks([-4 -3 -2 -1 0 1 2 3 4]);
    xticklabels(string([1E-4 1E-3 1E-2 1E-1 1 1E1 1E2 1E3 1E4]));
    title(tt);
    xline(0, 'LineStyle', '--', 'LineWidth', 2); grid on;

    % plot CDF
    [counts, edges] = histcounts(log10(valid), 'Normalization', 'pdf');

    % Bin widths
    bin_width = diff(edges);
    
    % CDF (integral of PDF)
    cdf_vals = cumsum(counts .* bin_width);
    
    % x-axis (bin centers)
    x_vals = edges(1:end-1) + bin_width/2;
    
    % Plot
    figure()
    plot(x_vals, cdf_vals, 'LineWidth', 1.5)
    grid on
    xlabel('relative error')
    ylabel('CDF '); title(tt);
    xticks([-4 -3 -2 -1 0 1 2 3 4]);
    xticklabels(string([1E-4 1E-3 1E-2 1E-1 1 1E1 1E2 1E3 1E4]));
end

