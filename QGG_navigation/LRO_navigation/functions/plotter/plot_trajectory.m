function [] = plot_trajectory(time,state)
    % Plot the spacecraft trajectory.
    
    r = state(:, 1:3)'; v = state(:, 4:6)';
    
    % Convert to datetime
    utc = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

    % plot orbit radius & 3D trajectory
    figure()
    subplot(1, 2, 1)
    plot(tUTC, vecnorm(r)./1E3, 'LineWidth', 2)
    grid on;
    title('orbit radius norm')
    xlabel('Epoch')
    ylabel('[km]')

    subplot(1, 2, 2)
    plot(tUTC, vecnorm(v), 'LineWidth', 2)
    grid on;
    title('orbit velocity norm')
    xlabel('Epoch')
    ylabel('[m/s]')
    hold all;

    figure()
    plot3(r(1, :), r(2, :), r(3, :), 'LineWidth', 2)
    hold on;
    
    % Create a sphere
    [R]  = cspice_bodvrd('MOON', 'RADII', 3)*1E3;
    R_M  = R(1);
    [Xs, Ys, Zs] = sphere(100);      % resolution 100x100
    surf(R_M*Xs, R_M*Ys, R_M*Zs, ...
        'FaceAlpha', 1, ...          % transparency
        'EdgeColor', 'none', ...
        'FaceColor', [0.8 0.8 0.8]); % light gray
    
    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('LRO Trajectory around the Moon. J2000 frame');

    axis equal

    % Plot orbit altitude
    figure();
    H = vecnorm(r) - R_M;   % [m]
    plot(tUTC, H./1E3, 'LineWidth', 2)
    grid on; xlabel('Epoch'); ylabel('[Km]');
    title('orbit altitude')
end

