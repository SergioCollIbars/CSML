function [Q_B] = rotate_processNoise(BN_mat,Q, k)
    %% ROTATE STATE

    % Initial time
    maxInd         = 3 * k; minInd = maxInd - 2;
    BN             = BN_mat(minInd:maxInd, :);

    A  = blkdiag(BN, BN, eye(6));
    
    Q_B = A * Q * A.';
end

