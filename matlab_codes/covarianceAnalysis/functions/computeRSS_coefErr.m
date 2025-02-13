function [CoefErr] = computeRSS_coefErr(n_max, Nc, Ns, err)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = (err(1))^2;
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefErr(n) = CoefErr(n)  + (err(j))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        CoefErr(n) = CoefErr(n)  + (err(j))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
    for n = 1:n_max
        CoefErr(n) = sqrt(CoefErr(n));
    end
end