function [X_B] = rotate_state(time,BN_mat,X_N)
    %% ROTATE STATE
    
    % Number of time steps
    Nt = length(time);

    X_B = X_N.*0;
    for k = 1:Nt
        maxInd         = 3 * k; minInd = maxInd - 2;
        BN             = BN_mat(minInd:maxInd, :);

        A = blkdiag(BN, BN);

        X_B(:, k)      = A * X_N(:, k);
    end
end

