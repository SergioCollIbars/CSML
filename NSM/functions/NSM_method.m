function [ax, nx] = NSM_method(Y, Yc, Hc, R, Hpr)
    dY = Y - Yc;
    C = null(Hpr');
    hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end); Hc(2, 2:end); Hc(5, 2:end);...
         Hc(8, 2:end);  Hc(3, 2:end); Hc(6, 2:end); Hc(9, 2:end)];

    % project measurements
    y  = C' * dY;
    hc = C' * hc;
    r  = C' * R * C;

    % information and normal matrices
    ax = hc' * (r\hc);
    nx = hc' * (r\y);
end

