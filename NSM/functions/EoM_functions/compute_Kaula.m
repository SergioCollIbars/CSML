function [sigmaK] = compute_Kaula(n_max, K)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigmaK = ones(1, Ncs);

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        sigmaK(j) = K/(n^2);
        
        % increment degree
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
        % compute degree variance
         sigmaK(j) = K/(n^2);

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end
