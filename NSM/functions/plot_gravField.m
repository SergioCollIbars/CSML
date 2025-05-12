function [] = plot_gravField(X, sigma_N, SH_N, n_max, tt, mk, llg)
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    figure()
    subplot(1, 2, 1)
    semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Nc-1, 3*sigma_N(1:Nc -1), 'Marker',mk, 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N(1:Nc-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    title(tt + "C_{nm}")
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    
    subplot(1, 2, 2)
    semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold on;
    semilogy(1:Ns, 3*sigma_N(Nc:Nc+Ns-1), 'Marker',mk, 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
    title(tt + "S_{nm}")
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    legend(llg)
end
