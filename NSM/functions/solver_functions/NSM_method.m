function [ax, nx, mxc, mcc] = NSM_method(Y, Yc, Hc, R, Hpr)
    dY = Y - Yc;
    % % C = null(Hpr');
    [~,~,d] = svd(Hpr');
    C = d(:, 6);

    % project measurements
    y  = C' * dY;
    hc = C' * Hc;
    r  = C' * R * C;
    hrot = C' * Hpr;
    
    % information and normal matrices
    ax = hc' * (r\hc);
    nx = hc' * (r\y);

    mxc = (hc' * (r\hrot));
    mcc = (hrot' * (r\hrot)); 
end

