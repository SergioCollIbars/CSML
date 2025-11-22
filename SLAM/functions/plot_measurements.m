function [] = plot_measurements(time,Y)
    % Description: plot measurements over time
    
    % Number of measurements
    Nm = length(Y(:, 1));
    
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    figure()
    for k = 1:Nm
        subplot(2, 3, k)
        plot(time./3600, Y(k, :)./1E3, 'LineWidth', 2)
        xlabel('Time [hr]')
        ylabel('[Eotvos]')
        grid on;
        title(tt(k));
    end
end

