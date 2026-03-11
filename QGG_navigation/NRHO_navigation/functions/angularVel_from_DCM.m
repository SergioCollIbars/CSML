function [angVel_vec] = angularVel_from_DCM(BN_mat, dt)
    % Givent the BN rotation matrix sequence with dimensions: 3 * Nt x 3,
    % obtain the angular velocity vector.
    % dt: time interval in seconds.

    Nt         = length(BN_mat(:, 1)) / 3;
    angVel_vec = zeros(3, Nt);
    for k = 2:Nt-1
        % current DCM
        maxInd = 3 * k; minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);
        
        % time series DCM
        i = k - 1;
        maxInd = 3 * i; minInd = maxInd - 2;
        BN_prev = BN_mat(minInd:maxInd, :);

        j = k + 1;
        maxInd = 3 * j; minInd = maxInd - 2;
        BN_next = BN_mat(minInd:maxInd, :);

        BN_dot = (BN_next - BN_prev)./(2*dt);
        
        % compute angular velocity
        omega_dyad = - BN_dot * BN';
        angVel_vec(:, k) = [omega_dyad(3, 2);omega_dyad(1,3);...
            omega_dyad(2, 1)];
    end
end

