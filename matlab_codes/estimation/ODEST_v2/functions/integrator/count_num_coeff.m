function [Nc, Ns, Ncs] = count_num_coeff(degree)
    % DESCRIPTION: count the number of zonal and sectorial / tesseral
    % coeffcients for a degree (n) value.
    Nc = 1;
    for k = 2:degree
        Nc = Nc + k + 1;
    end
    Ns = 0;
    for k = 2:degree
        Ns = Ns + k;
    end
    
    % total number of coefficients
    Ncs = Nc + Ns;
end

