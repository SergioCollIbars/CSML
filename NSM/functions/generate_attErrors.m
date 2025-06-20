function [att_Err] = generate_attErrors(t, type, planet, Amp, T, measMode)
    Nt  =length(t); 
    dt = t(2) - t(1);
    % define att_err vector
    att_Err = ones(6, Nt) * NaN;    % [att_err, datt_err]
    if(measMode == "2") % GG + SST X = (q, q_dot)
        if(type == "constant")
            att_Err(1:3, :)   = ones(3, Nt).*Amp;
            att_Err(4:6, :)   = zeros(3, Nt);
        elseif(type == "periodic")
            w = 1/2 * 2 * pi / T;
            att_Err(1:3, :) = Amp.*sin(w.*t);
            att_Err(4:6, :) = w*Amp.*cos(w.*t);
        elseif(type == "linear")
            time = t - t(1);
            periodic_time = mod(time, 3*T);  % Time within each orbital period
            att_Err(1:3, :) = Amp.*periodic_time;
            att_Err(4:6, :) = Amp.*ones(3, length(t));
        elseif(type == "random")
            att_Err(1:3, :) = normrnd(0, Amp(1), [3, length(t)]);
            att_Err(4:6, :) = centerDiff(att_Err(1:3, :), t(2)-t(1), Nt);
        end
    elseif(measMode == "1") % GG + SST + GYRO X = (q, omega)
        if(type == "constant")
            att_Err(1:3, :)   = ones(3, Nt).*Amp;
            att_Err(4:end, :) = ones(3, Nt).*Amp/dt;
        elseif(type == "periodic")
            w = 1/2 * 2 * pi / T;
            att_Err(1:3, :) = Amp.*sin(w.*t);
            att_Err(4:6, :) = w.*Amp.*cos(w.*t);
        elseif(type == "linear")
            time = t - t(1);
            periodic_time = mod(time, 3*T);  % Time within each orbital period
            att_Err(1:3, :) = Amp.*periodic_time;
            att_Err(4:6, :) = Amp/dt.*periodic_time;
        elseif(type == "random")
            att_Err(1:3, :) = normrnd(0, Amp(1), [3, length(t)]);
            % % att_Err(4:6, :) = normrnd(0, Amp(1)/dt, [3, length(t)]);
            att_Err(4:6, :) = centerDiff(att_Err(1:3, :), t(2)-t(1), Nt);
        end

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
    subplot(1, 2, 1)
    plot(t, att_Err(1:3, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad]')
    title('$\delta \theta$', 'Interpreter', 'latex')
    legend('\delta \Psi, yaw', '\delta \theta, pitch', '\delta \phi roll')
    grid on;

    subplot(1, 2, 2)
    plot(t, att_Err(4:6, :), 'LineWidth', 2)
    xlabel('Time')
    ylabel('[rad/s]')
    title('$\delta \dot{\theta}$', 'Interpreter', 'latex')
    grid on;

% %     subplot(1, 3, 3)
% %     plot(t, att_Err(7:9, :), 'LineWidth', 2)
% %     xlabel('Time')
% %     ylabel('[rad/s^2]')
% %     title('$\delta \ddot{\theta}$', 'Interpreter', 'latex')
% %     sgtitle('Attitude error');
% %     grid on;
end

function [df] = centerDiff(f, dt, Nt)
    df = f * NaN;
    for j = 2:Nt-1
        df(:, j) = (f(:, j+1) - f(:, j-1))./(2*dt);
    end
end
