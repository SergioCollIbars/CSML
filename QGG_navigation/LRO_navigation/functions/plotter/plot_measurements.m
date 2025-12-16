function [] = plot_measurements(time,Y, bias)
    % Description: plot measurements over time
    
    % Number of measurements
    Nm = length(Y(:, 1));
    
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    utc = cspice_et2utc(time', 'ISOC', 6);
    tUTC = datetime(utc, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

    figure()
    for k = 1:Nm
        subplot(2, 3, k)
        plot(tUTC, Y(k, :)./1E3, 'LineWidth', 2)
        xlabel('Epoch')
        ylabel('[Eotvos]')
        grid on;
        title(tt(k));
    end
    sgtitle('Gradiometer measurements');


    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    figure()
    for k = 1:Nm
        subplot(2, 3, k)
        plot(tUTC, bias(k, :), 'LineWidth', 2)
        xlabel('Epoch')
        ylabel('[milli-Eotvos]')
        grid on;
        title(tt(k));
    end
    sgtitle('Gradiometer bias');
end

