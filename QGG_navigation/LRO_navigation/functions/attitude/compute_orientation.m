function [BN_mat, NB_EARTH_mat, NB_MOON_mat] = ...
    compute_orientation(time, state, type)
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
    
    NB_EARTH_mat = ones(3*Nt, 3);
    NB_MOON_mat  = ones(3*Nt, 3);
    for k = 1:Nt
        et = time(k);
        % compute Earth and Moon body frame to J2000
        frame_to   = 'J2000';
        frame_from = 'IAU_EARTH';
        J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

        maxInd = 3 * k; minInd = maxInd - 2;
        NB_EARTH_mat(minInd:maxInd, :) = J2000_EARTH;
    
        frame_from = 'MOON_PA';
        J2000_MOON = cspice_pxform(frame_from, frame_to, et);
        NB_MOON_mat(minInd:maxInd, :)  = J2000_MOON;
    end
end

