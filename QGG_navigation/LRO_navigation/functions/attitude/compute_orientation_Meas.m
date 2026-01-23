function [Y_B] = compute_orientation_Meas(time, BN_mat, Y, signal_error, SF)
    % Description: compute the S/C orientation over time.

    % time vector length
    Nt = length(time);
    
    Y_B = Y.*0;
    for k =1:Nt
        maxInd = 3 * k; minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);

        % measurement inertial frame [mE]
        T_ACI = [Y(1, k),Y(2, k), Y(3, k);...
                 Y(2, k),Y(4, k), Y(5, k);...
                 Y(3, k),Y(5, k), Y(6, k)];
        T_B   = BN * T_ACI * BN';
    
        Y_B(:, k) = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);...
                     T_B(2,3);T_B(3,3)];
    end

    % Add error to true signal [mE]
    Y_B = SF.*Y_B + signal_error;
end

