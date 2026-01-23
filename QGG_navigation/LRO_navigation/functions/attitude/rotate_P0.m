function [P_B] = rotate_P0(BN_mat,P0)
    %% ROTATE STATE

    % Initial time
    k = 1;
    maxInd         = 3 * k; minInd = maxInd - 2;
    BN             = BN_mat(minInd:maxInd, :);

    A  = blkdiag(BN, BN, eye(6), eye(6));
    
    P_B = A * P0 * A.';
end

