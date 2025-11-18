function plot_measurements(t, T, planetParams, r_sc_moon)
    %%                    PLOT RESULTS FUNCTION
    % Description: Plot the measurements.
    % Author: Sergio Coll Ibars
    % Date: 03/29/2024
    
    % get indexes for local minima (prilune)
    idx_min = islocalmin(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations
    idx_max = islocalmax(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations


    % plot options
    lw = 3;
    color1 = [204, 0, 204]./256;     % violet
    color2 = 'g';
    set(0,'defaultAxesFontSize',16);

    % plot gradiometer measurements
    lbl = ["\Gamma_{xx} [E]", "\Gamma_{xy} [E]", "\Gamma_{xz} [E]", "\Gamma_{yy} [E]", ...
        "\Gamma_{yz} [E]", "\Gamma_{zz} [E]"];

    % select and format time
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

    figure()
    for j = 1:6
        subplot(3, 2, j)
        plot(time, T(j, :)*1E3, 'LineWidth', lw, ...
            'Color', color2, 'Marker','.')
        xlabel(xlb)
        ylabel(lbl(j))
        grid on;
        if(sum(idx_min)~=0)
            xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
                'LineStyle', '--')
        end
        if(sum(idx_max)~=0)
            xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
                'LineStyle', '--')
        end
    end
    sgtitle(tt)

    figure()
    plot(time, T(7:end, :)./1E3, ...
        'LineWidth', lw, 'Marker', '.')
    title('Accelerometer measurements')
    xlabel(xlb)
    ylabel('[m/s^2]')
    legend('a_x', 'a_y', 'a_z');
    grid on;

    T_norm = zeros(1, length(T(1, :)));
    for j = 1:length(T(1, :))
        TT = [T(1, j), T(2, j), T(3, j);...
            T(2, j), T(4, j), T(5, j);...
            T(3, j), T(5, j), T(6, j)];
        T_norm(j) = norm(TT, 'fro'); 
    end
    figure()
    plot(time, T_norm, 'LineWidth', 2, 'Color', color2);
    hold all;
    if(sum(idx_min)~=0)
        xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
        'LineStyle', '--')
    end
    if(sum(idx_max)~=0)
        xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
        'LineStyle', '--')
    end
    grid on;
    title('Gradiometer tensor norm');
end

