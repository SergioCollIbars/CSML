function [Ar] = generate_posErrors(t, type, Amp, T)
    Nt  =length(t);
    if(type == "constant")
        Ar = ones(3, Nt).*Amp;
    elseif(type == "periodic")
        w = 2 * pi / T;
        Ar = Amp.*sin(w.*t);
    elseif(type == "random")
        Ar = normrnd(0, Amp, [3, Nt]);
    end

    % plot position error
    figure()
    plot(t, Ar, 'LineWidth', 2)
    xlabel('Time')
    ylabel('[m]')
    title('position error. Inertial frame')
    legend('x', 'y', 'z')

end

