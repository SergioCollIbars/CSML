function [ax, nx, mxc, mcc] = LS_method(Y, Yc, hc, Hrot, R)
    dY = Y - Yc;
    hrot = Hrot;

    % information and normal matrices
    ax = hc' * (R\hc);
    nx = hc' * (R\dY);
    mxc = (hc' * (R\hrot));
    mcc = (hrot' * (R\hrot)); 
end