function [] = do_smoothing(time, state_t, planetParams, poleParams,...
    instrumentParams, BN_mat, Cnm_t, Snm_t, Y, Xs, Ps, Xg, Pg)

     % load filter parameters
    [R0, ~,  ~, Q0, Qb, ~, ~, ~, smoothing] = ...
        loadFilterParams(planetParams, instrumentParams, Cnm_t, Snm_t);

    if(smoothing)
        disp('Applying smoothing');
        Ns = length(Xs(:, end));

        % assign initial values (start at final time)
        P0 = reshape(Ps(end, :), [Ns, Ns]);
        P0c_grav = Pg;
        
        [Ncoeff, Nscoeff, Ncs]  = count_num_coeff(planetParams(3)); 
        [Cnm, Snm] = list2mat(planetParams(3), Ncoeff, Nscoeff, Xg);

        % run EKF
        [Xs_smooth, Ps_smooth] = EKF_smoothing(time, planetParams, ...
            poleParams, Y, Xs(:, end), P0, ...
            P0c_grav, BN_mat, Cnm, Snm, Q0, Qb, R0);
        
        [std_pos, std_bias] = plot_states(Xs, Ps, Xs_smooth, ...
            Ps_smooth, time, state_t);
        
        % update gravity estimation
        Y_cal     = Y - Xs_smooth(7:end, :);
        sigmaPos  = max(max(std_pos(:, 1:end-1)'));
        sigmaBias = max(std_bias(:, 1:end-1)');

        [Xg_smooth, Pg_smooth] = gravEstim_smoothing(time, Y_cal, planetParams,...
            poleParams, Xs_smooth, P0c_grav, ...
            Cnm, Snm, sigmaPos, sigmaBias, R0);
        
        Xc_t   = mat2list(Cnm_t, Snm_t, Ncoeff, Nscoeff);
        err_smth = abs(Xc_t - Xg_smooth);
        err      = abs(Xc_t - Xg);

        std_smth = sqrt(diag(Pg_smooth));
        std      = sqrt(diag(Pg));

        figure()
        semilogy(1:Ncs-1, abs(Xc_t(2:end)), 'LineWidth', 2, 'Color', 'b');
        hold all;
        semilogy(1:Ncs-1, err(2:end), 'LineStyle', 'none', 'Color', ...
            'r', 'Marker', '*');
        semilogy(1:Ncs-1, 3.*std(2:end), 'LineWidth', 2, 'Color', 'k');

        semilogy(1:Ncs-1, err_smth(2:end), 'LineStyle', 'none', 'Color', ...
            'g', 'Marker', '*');
        semilogy(1:Ncs-1, 3.*std_smth(1:end), 'LineWidth', 2, 'Color', 'g');

% %         Xg = Xg_smooth; 
% %         Pg = zeros(Ncs);
% %         Pg(2:end, 2:end) = Pg_smooth;
    end

    % save gravity field uncertainty
    disp('Save Gravity data ...');
    save("data/gravField_cov.mat", 'Pg');
    save("data/gravField_coeff.mat", 'Xg');
    
    
end


%% FUNCTIONS
function [std_pos, std_bias] = plot_states(Xs, Ps, Xs_smooth,...
    Ps_smooth, time, state_t)
    Nt = length(time); Ns = length(Xs(:, 1));

    std = nan(Ns, Nt); std_smth = nan(Ns, Nt);
    for k = 1:Nt
        p = reshape(Ps(k, :), [Ns, Ns]);
        std(:, k) = sqrt(diag(p)); 

        p_smth = reshape(Ps_smooth(k, :), [Ns, Ns]);
        std_smth(:, k) = sqrt(diag(p_smth)); 
    end

    err      = abs(Xs - state_t(:, 1:Ns)');
    err_smth = abs(Xs_smooth - state_t(:, 1:Ns)');
    
    % plot position norm error
    figure()
    subplot(1,2,1)
    semilogy(time./3600, vecnorm(err(1:3, :)), 'Color', 'r'); hold all;
    semilogy(time./3600, vecnorm(err_smth(1:3, :)), 'Color', 'g');
    semilogy(time./3600, 3.*vecnorm(std(1:3, :)), 'Color', 'k', 'LineWidth', 2);
    semilogy(time./3600, 3.*vecnorm(std_smth(1:3, :)), 'Color', 'g', 'LineWidth', 2);
    grid on;
    ylabel('[m]'); xlabel('[hr]');
    title('position');

    subplot(1,2,2)
    semilogy(time./3600, vecnorm(err(4:6, :)), 'Color', 'r'); hold all;
    semilogy(time./3600, vecnorm(err_smth(4:6, :)), 'Color', 'g');
    semilogy(time./3600, 3.*vecnorm(std(4:6, :)), 'Color', 'k', 'LineWidth', 2);
    semilogy(time./3600, 3.*vecnorm(std_smth(4:6, :)), 'Color', 'g', 'LineWidth', 2);
    grid on;
    ylabel('[m]'); xlabel('[hr]');
    title('velocity')
    sgtitle('EKF smoothing results');

    % retur STD values for position and bias
    std_pos  = std_smth(1:3, :);
    std_bias = std_smth(7:end, :);
end
