function [C_mat,S_mat] = get_SHmatrices(estimObj, X)
     % compute Cmat and Smat at nominal
    C_mat = zeros(estimObj.n_max + 1, estimObj.n_max + 1);
    S_mat = zeros(estimObj.n_max + 1, estimObj.n_max + 1);
    Nxc = estimObj.Nc;
    Nxs = estimObj.Ns;
    
    C_mat(1, 1) = 1;
    
    n = 2;
    m = 0;
    for j = 2:Nxc
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
    for j = Nxc + 1:Nxs + Nxc
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

