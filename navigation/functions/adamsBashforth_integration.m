function [x2] = adamsBashforth_integration(x, f, At, order)
    %%              ADAMS BASHFORTH INTEGRATION METHOD
    % Description: integration based on two steps Adams bashforth. Returns the
    % variable at x(2) given X(0), X(1), its time derivative f(0), f(1) and the time
    % step At
    
    if(order == 2)
        x1 = reshape(x, [6, 6]);
        f1 = reshape(f(1:36), [6, 6]);
        f0 = reshape(f(37:end), [6, 6]);
        x2  = x1 + 3/2 * At * f1 - 0.5 * At * f0;
    elseif(order == 3)
        x1 = reshape(x, [6, 6]);
        f2 = reshape(f(1:36), [6, 6]);
        f1 = reshape(f(37:(2*36)), [6, 6]);
        f0 = reshape(f((2*36+1):end), [6, 6]);
        x2 = x1 + At * (23/12*f2 - 16/12*f1 + 5/12*f0);
    elseif(order == 4)
        x1 = reshape(x, [6, 6]);
        f3 = reshape(f(1:36), [6, 6]);
        f2 = reshape(f(37:(2*36)), [6, 6]);
        f1 = reshape(f((2*36+1):36*3), [6, 6]);
        f0 = reshape(f((3*36+1):end), [6, 6]);
        x2 = x1 + At * (55/24*f3 - 59/24*f2 + 37/24*f1 - 9/24*f0);
    elseif(order == 5)
        x1 = reshape(x, [6, 6]);
        f4 = reshape(f(1:36), [6, 6]);
        f3 = reshape(f(37:(2*36)), [6, 6]);
        f2 = reshape(f((2*36+1):36*3), [6, 6]);
        f1 = reshape(f((3*36+1):4*36), [6, 6]);
        f0 = reshape(f((4*36+1):end), [6, 6]);
        x2 = x1 + At * (1901/720*f4 - 2774/720*f3 + 2616/720*f2 - ...
            1274/720*f1 + 251/720*f0);
end

