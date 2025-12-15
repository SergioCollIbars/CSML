function [] = plot_trajectory(time,state)
    % Plot the spacecraft trajectory.
    
    r = state(:, 1:3)'; v = state(:, 4:6)';
    
    % plot orbit radius & 3D trajectory
    figure()
    subplot(1, 2, 1)
    plot(time./3600, vecnorm(r)./1E3, 'LineWidth', 2)
    grid on;
    title('orbit radius norm')
    xlabel('Time [hr]')
    ylabel('[km]')

    subplot(1, 2, 2)
    plot3(r(1, :)./1E3, r(2, :)./1E3, r(3, :)./1E3, 'LineWidth', 2)
    grid on;
    title('3D trajectory ACI frame')
    xlabel('X [Km]')
    ylabel('Y [km]')
    zlabel('Z [Km]')
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
    
    axis equal;
    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    title('LRO Trajectory around the Moon. J2000 frame')

    axis equal
end

