function [Xg_smth, Pg_smth, Xs_smooth, Ps_smooth] = do_smoothing(time, state_t, planetParams,...
    poleParams,instrumentParams, BN_mat, Cnm_t, Snm_t, Y, Xs, Ps,...
    Xg, Pg, h, ax)

     % load filter parameters
    [R0, ~,  P0c, Q0, Qb, ~, ~, ~, ~, smoothing] = ...
        loadFilterParams(planetParams, instrumentParams, Cnm_t, Snm_t);

    if(smoothing == 1)
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
            (1.3).*P0c_grav, BN_mat, Cnm, Snm, Q0, Qb, R0);
        Xs_smooth(:, end) = Xs(:, end);
        
        % compute smooth errors
        err_smooth = abs(Xs_smooth - state_t(:, 1:12)');

% %         % plot smoothing results
% %         [h, ax] = plot_smooth_results(time, h, ax, ...
% %             err_smooth, Ps_smooth);

        % update gravity estimation
        Y_cal     = Y - Xs_smooth(7:end, :);

        [Xg_smooth, Pg_smooth] = gravEstim_smoothing(time, Y_cal, ...
            planetParams, poleParams, Xs_smooth, P0c, ...
            Cnm, Snm, R0, Ps_smooth);
        
        Xc_t   = mat2list(Cnm_t, Snm_t, Ncoeff, Nscoeff);
        errC_smth = abs(Xc_t - Xg_smooth);

        stdC_smth = sqrt(diag(Pg_smooth));

% %         green = [0 1 0];   % MATLAB default green
% %         dark_green = green * 0.5;
% % 
% %         [errC_order, xvals] = orderValues(errC_smth(2:end), planetParams(3));
% %         [sigmaC_order, ~]   = orderValues(stdC_smth, planetParams(3));
% % 
% %         % update gravity estimation figure 
% %         gray = 0.7;
% %         set(h.Err_c, 'Color', [gray, gray, gray]);
% %         gray = 0.3;
% %         set(h.Cov_c, 'Color', [gray, gray, gray]);
% % 
% %         semilogy(ax.grav, xvals, abs(errC_order),'LineStyle', 'none', ...
% %             'Marker', '*', "MarkerSize", 2, 'Color', green)
% %         semilogy(ax.grav, xvals, 3.*sigmaC_order,'LineStyle', '-', ...
% %             'LineWidth', 2, "MarkerSize", 2, 'Color', dark_green)
    else
        [~, ~, Ncs]  = count_num_coeff(planetParams(3)); 
        Pg_smooth    =  Pg(2:end, 2:end);

        Xg_smooth    =  Xg;
    end

    Pg_smth = zeros(Ncs);
    Pg_smth(2:end, 2:end) = Pg_smooth;

    Xg_smth = Xg_smooth; 
end
