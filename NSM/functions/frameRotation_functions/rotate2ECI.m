function [state_ACI] = rotate2ECI(state_ECEF, ACI_ECEF, t)
    state_ACI = ones(6, length(t)) * NaN;
    Nt  = length(t);
    R_dot = ones(3,3,Nt) * NaN;
    for j = 2:Nt-1
        dt = t(j+1) - t(j-1);
        
        maxPos = 3*(j-1); minPos = maxPos- 2;
        R1 = ACI_ECEF(minPos:maxPos, :);

        maxPos = 3*(j+1); minPos = maxPos- 2;
        R2 = ACI_ECEF(minPos:maxPos, :);

        R_dot(:,:, j) = (R2 - R1)./dt;
    end

    for j = 1:length(t)
        % rotation matrix
        maxPos = 3*j; minPos = maxPos- 2;
        R = ACI_ECEF(minPos:maxPos, :);
        Rd = R_dot(:, :, j);

        state_ACI(1:3, j) = R * state_ECEF(1:3, j);
        state_ACI(4:6, j) = Rd * state_ECEF(1:3, j) + R * state_ECEF(4:6, j);
    end

end