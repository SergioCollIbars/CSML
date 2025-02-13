function [C_mat, S_mat] = list2mat(n_max, Nc, Ns, X)
    
    % define output matrices
    C_mat = zeros(n_max + 1, n_max + 1);
    S_mat = zeros(n_max + 1, n_max + 1);
    
    % fill matrices
    C_mat(1, 1) = 1;
    
    n = 2;
    m = 0;
    for j = 2:Nc
        N = n + 1;
        M = m + 1;
        C_mat(N, M) = X(j);
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n +1;
        end
    end
    
    n = 2;
    m = 0;
    for j = Nc + 1:Ns + Nc
        N = n + 1;
        M = m + 2;
        S_mat(N, M) = X(j);
        if(m < n - 1)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end
end

