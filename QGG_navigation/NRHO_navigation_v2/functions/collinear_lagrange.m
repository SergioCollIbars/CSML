function y = collinear_lagrange(xstar)
        m1 = 5.974E24; % kg
        m2 = 7.348E22; % kg
        pi2 = m2/(m1 + m2);
        firstterm = xstar;
        secondterm = (1 - pi2) ./ abs(xstar + pi2).^3 .* (xstar + pi2);
        thirdterm = pi2 ./ abs(xstar - 1 + pi2).^3 .* (xstar - 1 + pi2);
        y = firstterm - secondterm - thirdterm;
end
