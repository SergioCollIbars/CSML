function [N] = NormFactor(n, m)
    % Description: given degree, n and order, m compute the normalice
    % factor
    if(m == 0)
        delta = 1;
    else
        delta = 0;
    end
    fac1 = factorial(n - m);
    fac2 = factorial(n + m);
    N = ((2 - delta)*(2*n + 1) * fac1 /fac2)^(0.5);
end
