function [errNorm_pos, errNorm_vel, ...
    errNorm_bias, std_pos, std_vel, std_bias] = ...
    update_plots(time, k, state_true, X0, P, h, ax, errNorm_pos,...
    errNorm_vel, errNorm_bias, std_pos, std_vel, std_bias, errC, sigmaC, Ncs)
        % state and bias errors
        e_k            = state_true(k, 1:6)'  - X0(1:6);
        e_b            = state_true(k, 7:12)' - X0(7:12); 
        errNorm_pos(k) = norm(e_k(1:3)); errNorm_vel(k) = norm(e_k(4:6));
        errNorm_bias(:, k) = abs(e_b);

        % state and bias std
        std = diag(sqrt(P));
        std_pos(k)     = 3*norm(std(1:3)); std_vel(k)   = 3*norm(std(4:6));
        std_bias(:, k) = 3.*std(7:12);

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

        set(h.Err_c, 'XData', 1:Ncs-1, 'YData',...
        errC(2:end), 'LineWidth', 1.5, 'Color', 'r');
        set(h.Cov_c, 'XData', 1:Ncs-1, 'YData', 3.*sigmaC(2:end), ...
        'LineWidth', 2, 'Color', 'k');
    
        drawnow limitrate;  % update the figure without killing performance
end

