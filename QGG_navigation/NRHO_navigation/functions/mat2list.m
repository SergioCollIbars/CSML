function [X] = mat2list(C_mat, S_mat, Nxc, Nxs)
    X = zeros(Nxc + Nxs, 1);
    X(1) = 1;
    
    n = 2;
    m = 0;
    for j = 2:Nxc
        N = n + 1;
        M = m + 1;
        X(j) = C_mat(N, M);
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
        X(j) = S_mat(N, M);
        if(m < n - 1)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end  
end

