function [ax, nx, mxc, mcc] = LS_method(Y, Yc, Hc, Hrot, R)
    dY = Y - Yc;
    hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end); Hc(2, 2:end); Hc(5, 2:end);...
             Hc(8, 2:end);  Hc(3, 2:end); Hc(6, 2:end); Hc(9, 2:end)];
    hrot = Hrot;

    % information and normal matrices
    ax = hc' * inv(R) * hc;
    nx = hc' * inv(R) * dY;
    mxc = (hc' * inv(R) * hrot);
    mcc = (hrot' * inv(R) * hrot); 
end