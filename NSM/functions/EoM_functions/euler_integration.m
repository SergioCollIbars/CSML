function [x] = euler_integration(x0, f, At)
    %%              EULER INTEGRATION METHOD
    % Description: integration based on finite differences. Returns the
    % variable at x(t) given X(0), its time derivative f(0) and the time
    % step At

    x = x0 + f * At;
end

