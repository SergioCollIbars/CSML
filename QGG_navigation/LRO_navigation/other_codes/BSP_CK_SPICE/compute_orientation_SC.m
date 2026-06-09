function [BN_mat] = compute_orientation_SC(time, state, type, NB_MOON)
    % Description: compute the S/C orientation over time.

    % time vector length
    Nt = length(time);

    BN_mat = ones(3*Nt, 3);
    if(type == "Inertial")
        for k =1:Nt
            maxInd = 3 * k; minInd = maxInd - 2;
            BN_mat(minInd:maxInd, :) = eye(3);
        end
    elseif(type == "RTN")
        for k =1:Nt
            r = state(k, 1:3)'; v = state(k, 4:6)';
            NB = RTN2ECI(r, v); % from RNT to ACI

            maxInd = 3 * k; minInd = maxInd - 2;
            BN_mat(minInd:maxInd, :) = NB';
        end
    elseif(type == "ENU")
        for k =1:Nt
            maxInd = 3 * k; minInd = maxInd - 2;
            J2000_ECEF = NB_MOON(minInd:maxInd, :);
            r_N        = state(k, 1:3)';
            X0_ECEF    = J2000_ECEF' * r_N;
            ENU_ECEF   = ecef2enu(X0_ECEF);
            BN         = ENU_ECEF * J2000_ECEF';

            BN_mat(minInd:maxInd, :) = BN;
        end
    else
        error("Orientation not defined.")
    end
end

