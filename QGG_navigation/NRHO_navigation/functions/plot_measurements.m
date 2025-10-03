function plot_measurements(t, T, planetParams, augmented_st, system)
    %%                    PLOT RESULTS FUNCTION
    % Description: Plot the measurements.
    % Author: Sergio Coll Ibars
    % Date: 03/29/2024
    
    % change tensor units to Etvos
    measDim_QGG =  (planetParams(3)^2*1E9);                 % [Etvos]
    measDim_Acc = (planetParams(2) * planetParams(3)^2);    % [m/s^2]

    % plot options
    lw = 3;
    color1 = [204, 0, 204]./256;     % violet
    color2 = 'g';
    set(0,'defaultAxesFontSize',16);

    % plot gradiometer measurements
    lbl = ["\Gamma_{11} [E]", "\Gamma_{12} [E]", "\Gamma_{13} [E]", "\Gamma_{22} [E]", ...
        "\Gamma_{23} [E]", "\Gamma_{33} [E]"];

    % select and format time
    if(system == "EPHEM")
        jd = 2451545 + t / planetParams(3) / 86400;
        humanReadableTime = datetime(jd, 'ConvertFrom', ...
            'juliandate');
        humanReadableTime.Format = 'MMM dd, yyyy';
        date_init = string(humanReadableTime(1));
        date_end  = string(humanReadableTime(end));
        humanReadableTime.Format = 'MMM dd';

        time = humanReadableTime;
        xlb = "date";
        tt = 'Gradiometer measurements from '  + date_init + ' - ' + date_end;
    else
        time = t/planetParams(3)/86400;
        xlb = "days";
        tt = 'Gradiometer measurements in CR3BP system';
    end

    figure()
    for j = 1:6
        subplot(3, 2, j)
        plot(time, T(j, :).*measDim_QGG, 'LineWidth', lw, ...
            'Color', color2, 'Marker','.')
        xlabel(xlb)
        ylabel(lbl(j))
        grid on;
    end
    sgtitle(tt)

    if(augmented_st)
        figure()
        plot(time, T(7:end, :).*measDim_Acc, ...
            'LineWidth', lw, 'Marker', '.')
        title('Accelerometer measurements')
        xlabel(xlb)
        ylabel('[m/s^2]')
        legend('a_x', 'a_y', 'a_z');
        grid on;
    end
end

