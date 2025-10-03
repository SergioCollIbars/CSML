function [Nc, Ns] = countCoeff(n)
    % count number of C and S coeff
    Nc = 1;
    for k = 2:n
        Nc = Nc + k + 1;
    end
    Ns = 0;
    for k = 2:n
        Ns = Ns + k;
    end
end

