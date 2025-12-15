function [h, ax] = plot_smooth_results(time, h, ax,...
    err_smooth, Ps_smooth)
    Nt = length(time); Ns = length(err_smooth(:, 1));

    std_smth = nan(Ns, Nt);
    for k = 1:Nt
        p_smth = reshape(Ps_smooth(k, :), [Ns, Ns]);
        std_smth(:, k) = sqrt(diag(p_smth)); 
    end

    % selec time frame
    th = time./(3600);

    green = [0 1 0];   % MATLAB default green
    dark_green = green * 0.5;

    % update position and velocity plots
    gray = 0.7;
    set(h.Err_pos, 'Color', [gray, gray, gray]);
    gray = 0.3;
    set(h.Cov_pos, 'Color', [gray, gray, gray]);
    semilogy(ax.pos, th, vecnorm(err_smooth(1:3, :)), 'LineWidth', 1.5, ...
        'Color', green);
    semilogy(ax.pos, th, 3.*vecnorm(std_smth(1:3, :)), 'LineWidth', 2, ...
        'Color', dark_green);

    gray = 0.7;
    set(h.Err_vel, 'Color', [gray, gray, gray]);
    gray = 0.3;
    set(h.Cov_vel, 'Color', [gray, gray, gray]);
    semilogy(ax.vel, th, vecnorm(err_smooth(4:6, :)).*(1E3), 'LineWidth', 1.5, ...
        'Color', green);
    semilogy(ax.vel, th, 3.*vecnorm(std_smth(4:6, :)).*(1E3), 'LineWidth',...
        2, 'Color', dark_green);

    for f = 1:6
            gray = 0.7;
            set(h.BiasErr(f), 'Color', [gray, gray, gray]);
            gray = 0.3;
            set(h.BiasCov(f), 'Color', [gray, gray, gray]);
            semilogy(ax.bias(f), th, err_smooth(6+f, :), 'LineWidth', 1.5, ...
            'Color', green);
            semilogy(ax.bias(f), th, 3.*std_smth(6+f, :), 'LineWidth', 2, ...
           'Color', dark_green);
    end
end

