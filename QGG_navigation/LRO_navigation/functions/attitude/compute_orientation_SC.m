function [BN_mat] = compute_orientation_SC(time, state, type)
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
    else
        error("Orientation not defined.")
    end
end

