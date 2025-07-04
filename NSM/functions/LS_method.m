function [ax, nx, mxc, mcc] = LS_method(Y, Yc, Hc, Hrot, R)
    dY = Y - Yc;
    hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end); Hc(2, 2:end); Hc(5, 2:end);...
             Hc(8, 2:end);  Hc(3, 2:end); Hc(6, 2:end); Hc(9, 2:end)];
    hrot = Hrot;

    % information and normal matrices
    ax = hc' * (R\hc);
    nx = hc' * (R\dY);
    mxc = (hc' * (R\hrot));
    mcc = (hrot' * (R\hrot)); 
end