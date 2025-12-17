function [] = plot_measurements(time,Y, bias, instrumentParams)
    % Description: plot measurements over time
    
    % measurement mask
    mask = instrumentParams(:, 1);
    gray = [.7 .7 .7];

    % Number of measurements
    Nm = length(Y(:, 1));
    
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    utc = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

    figure()
    for k = 1:Nm
        subplot(2, 3, k)
        if(mask(k)==1), color = 'b'; else, color = gray; end
        plot(tUTC, Y(k, :)./1E3, 'LineWidth', 2, 'Color', color);
        ylabel('[Eotvos]')
        grid on;
        title(tt(k));
    end
    sgtitle('Gradiometer measurements');


    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    figure()
    for k = 1:Nm
        subplot(2, 3, k)
        if(mask(k)==1), color = 'b'; else, color = gray; end
        plot(tUTC, bias(k, :), 'LineWidth', 2, 'Color', color)
        ylabel('[milli-Eotvos]')
        grid on;
        title(tt(k));
    end
    sgtitle('Gradiometer bias');
end

