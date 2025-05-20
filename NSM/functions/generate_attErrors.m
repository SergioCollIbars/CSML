function [att_Err] = generate_attErrors(t, type, planet, Amp, T)
    Nt  =length(t); 

    % define att_err vector
    att_Err = ones(9, Nt) * NaN;    % [att_err, datt_err, ddatt_err]
    if(type == "constant")
        att_Err(1:3, :)   = ones(3, Nt).*Amp;
        att_Err(4:end, :) = zeros(6, Nt);
    elseif(type == "periodic")
        w = 1/2 * 2 * pi / T;
        att_Err(1:3, :) = Amp.*sin(w.*t);
        att_Err(4:6, :) = w*Amp.*cos(w.*t);
        att_Err(7:9, :) = -w^2*Amp.*sin(w.*t);
    end

    % plot attitude error
    if(planet == "Earth")
         % convert GPS time to UTC
         gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
         t_UTC = gps_epoch + seconds(t);        % date time
         
         plot_atittude_err(t_UTC, att_Err);
         xlim_limits = [datetime(2012,11,02,0,0,0), datetime(2012,11,02,3,0,0)];
         set(findall(gcf, 'Type', 'axes'), 'XLim', xlim_limits)
    else
        plot_atittude_err(t, att_Err);
    end


end


function [] = plot_atittude_err(t, att_Err)
    figure()
    subplot(1, 3, 1)
    plot(t, att_Err(1:3, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad]')
    title('$\delta \theta$', 'Interpreter', 'latex')
    legend('\delta \Psi, yaw', '\delta \theta, pitch', '\delta \phi roll')
    grid on;

    subplot(1, 3, 2)
    plot(t, att_Err(4:6, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad/s]')
    title('$\delta \dot{\theta}$', 'Interpreter', 'latex')
    grid on;

    subplot(1, 3, 3)
    plot(t, att_Err(7:9, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad/s^2]')
    title('$\delta \ddot{\theta}$', 'Interpreter', 'latex')
    sgtitle('Attitude error');
    grid on;
end
