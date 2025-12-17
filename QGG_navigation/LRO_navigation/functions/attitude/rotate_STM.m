function [STM_B_vec] = rotate_STM(time,BN_mat,STM_vec, BN0)
    %% ROTATE STATE
    
    % Number of time steps
    Nt = length(time);
    Ns = sqrt(length(STM_vec(1, :)));
    
    STM_B_vec = STM_vec.*0;
    for k = 1:Nt
        maxInd         = 3 * k; minInd = maxInd - 2;
        BN             = BN_mat(minInd:maxInd, :);

        STM_N = reshape(STM_vec(k, :), [Ns, Ns]);
        A  = blkdiag(BN, BN);
        A0 = blkdiag(BN0, BN0); 

        STM_B = A * STM_N / A0;
        STM_B_vec(k, :) = reshape(STM_B, [1, Ns*Ns]);
    end
end

