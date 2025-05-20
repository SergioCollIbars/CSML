function [rn_ECI] = rotate2ECI(rn_ECEF, ACI_ECEF, t)
    rn_ECI = ones(3, length(t)) * NaN;

    for j = 1:length(t)
        % rotation matrix
        maxPos = 3*j; minPos = maxPos- 2;
        R = ACI_ECEF(minPos:maxPos, :);

        rn_ECI(:, j) = R * rn_ECEF(:, j);
    end
end