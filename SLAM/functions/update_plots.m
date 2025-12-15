function [errNorm_pos, errNorm_vel, ...
    errNorm_bias, std_pos, std_vel, std_bias] = ...
    update_plots(time, k, state_true, X, P, h, ax, errNorm_pos,...
    errNorm_vel, errNorm_bias, std_pos, std_vel, std_bias, errC, sigmaC, n_max)
        % state and bias errors
        e_k            = state_true(1:k, 1:6)'  - X(1:6, 1:k);
        e_b            = state_true(1:k, 7:12)' - X(7:12, 1:k); 
        errNorm_pos(1:k) = vecnorm(e_k(1:3, :)); 
        errNorm_vel(1:k) = vecnorm(e_k(4:6, :));
        errNorm_bias(:, 1:k) = abs(e_b);

        % state and bias std
        for j = 1:k
            std = diag(sqrt(reshape(P(j, :), [12, 12])));
            std_pos(j)     = 3*norm(std(1:3)); std_vel(j)   = 3*norm(std(4:6));
            std_bias(:, j) = 3.*std(7:12);
        end

        % select time scale
        th = time(1:k)./3600;
        maxScale = 4;

        % update X and Y data of the existing line objects
        set(h.Err_pos, 'XData', th, 'YData',...
            errNorm_pos(1:k), 'LineWidth', 1.5, 'Color', 'r');
        set(h.Cov_pos, 'XData', th, 'YData', std_pos(1:k), ...
            'LineWidth', 2, 'Color', 'k');
        set(ax.pos, 'YLim', [0, maxScale*std_pos(k)]);

        set(h.Err_vel, 'XData', th, 'YData',...
            errNorm_vel(1:k).*1E3, 'LineWidth', 1.5, 'Color', 'r');
        set(h.Cov_vel, 'XData', th, 'YData', std_vel(1:k).*1E3, ...
            'LineWidth', 2, 'Color', 'k');
        set(ax.vel, 'YLim', [0, maxScale*std_vel(k)*1E3]);

        for f = 1:6
            set(h.BiasErr(f), 'XData', th, 'YData',...
            errNorm_bias(f, 1:k), 'LineWidth', 1.5, 'Color', 'r');
            set(h.BiasCov(f), 'XData', th, 'YData', std_bias(f, 1:k), ...
            'LineWidth', 2, 'Color', 'k');
            set(ax.bias(f), 'YLim', [0, maxScale*std_bias(f, k)]);
        end
        
        % order values
        [errC_order, xvals] = orderValues(errC(2:end), n_max);
        [sigmaC_order, ~] = orderValues(sigmaC(2:end), n_max);

        set(h.Err_c, 'XData', xvals, 'YData',...
        errC_order, 'Color', 'r', 'Marker','*',...
            "MarkerSize", 2, 'LineStyle','none');
        set(h.Cov_c, 'XData', xvals, 'YData', 3.*sigmaC_order, ...
            'LineWidth', 2, 'Color', 'r');
    
        drawnow limitrate;  % update the figure without killing performance
end

