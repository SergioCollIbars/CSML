function [] = plot_SC_state(t, planet, state_t, orientation)
    % Description: plot the position and attitude orientation of the
    % spacecraft over time. 
    
    % plot for Earth
    if(planet == "Earth")
        % convert GPS time to UTC
        gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
        t_UTC = gps_epoch + seconds(t);        % date time 
    
        % plot position
        plot_position(t_UTC, state_t);
        xlim([datetime(2012,11,02,0,0,0), datetime(2012,11,02,3,0,0)])
        
        % plot velocity
        plot_velocity(t_UTC, state_t);
        xlim([datetime(2012,11,02,0,0,0), datetime(2012,11,02,3,0,0)])

        % plot attitude
        plot_EulerAngles(t_UTC, orientation)
        xlim_limits = [datetime(2012,11,02,0,0,0), datetime(2012,11,02,3,0,0)];
        set(findall(gcf, 'Type', 'axes'), 'XLim', xlim_limits)
    
    else
        % plot position
        plot_position(t, state_t);
       
        % plot velocity
        plot_velocity(t, state_t);

        % plot attitude
        plot_EulerAngles(t, orientation)
    end
end


% plot position
function [] = plot_position(t, state_t)
    figure()
    plot(t, state_t(:, 1:3)', 'LineWidth', 2)
    xlabel('Time')
    ylabel('[m]')
    legend('X', 'Y', 'Z')
    grid on;
    title('S/C position over time')
end

% plot velocity
function [] = plot_velocity(t, state_t)
    figure()
    plot(t, state_t(:, 4:6)', 'LineWidth', 2)
    xlabel('Time')
    ylabel('[m / s]')
    legend('X', 'Y', 'Z')
    grid on;
    title('S/C velocity over time')
end

% plot attitude (Euler angles)
function [] = plot_EulerAngles(t, orientation)
    figure()
    subplot(1, 3, 1)
    plot(t, orientation(1:3, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad]')
    grid on;
    legend({'$\Psi, yaw$', '$\theta, pitch$', '$\phi, roll$'}, ...
        'Interpreter', 'latex')

    subplot(1, 3, 2)
    plot(t, orientation(4:6, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad / s]')
    grid on;
    legend({'$\dot{\Psi}, yaw$', '$\dot{\theta}, pitch$', '$\dot{\phi}, roll$'}, ...
        'Interpreter', 'latex')

    subplot(1, 3, 3)
    plot(t, orientation(7:9, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad / s^2]')
    grid on;
    legend({'$\ddot{\Psi}, yaw$', '$\ddot{\theta}, pitch$', '$\ddot{\phi}, roll$'}, ...
        'Interpreter', 'latex')

    sgtitle('S/C orientation')
end
