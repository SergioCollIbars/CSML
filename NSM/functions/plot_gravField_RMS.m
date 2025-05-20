function [] = plot_gravField_RMS(X, sigma_N, SH_N, n_max, tt, llg)
    [Nc, Ns, ~] = count_num_coeff(n_max); 

    % compute RMS value true
     X_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

     % compute error RMS
     err_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                X - [1;SH_N], zeros(n_max+1, n_max+1), ...
                zeros(n_max+1, n_max+1)); 

     % compute uncertaity RMS (3 sigma)
     s_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
               [1;3*sigma_N], zeros(n_max+1, n_max+1), ...
               zeros(n_max+1, n_max+1));
    
     % plot RMS value
    figure()
    semilogy(1:n_max, X_RMS, 'LineWidth', 2, 'Color', 'k');
    hold all;
    semilogy(1:n_max, s_RMS, 'LineWidth', 2, 'Color', 'b');
    semilogy(1:n_max, err_RMS, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--');
    title(tt); legend(llg);
    xlim([2, n_max]);

end

