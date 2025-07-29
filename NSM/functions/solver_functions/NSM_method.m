function [ax, nx] = NSM_method(Y, Yc, Hc, R, Hpr)
    dY = Y - Yc;
    C = null(Hpr');

    % project measurements
    y  = C' * dY;
    hc = C' * Hc;
    r  = C' * R * C;
    
    % information and normal matrices
    ax = hc' * (r\hc);
    nx = hc' * (r\y);
end

