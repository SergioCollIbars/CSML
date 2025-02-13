function [ddU_ACI, H_ACI, H_RTN] = gradiometer_meas(time ,asterParams, poleParams, state, ...
                noise, Cnm, Snm)
    % GRADIOMETER_MEAS generate gradiometer measurement
    % Given the set of Cnm and Snm coefficient, the asteroid pole
    % parameters and the spacecraft position at different times, generate
    % gradiometer measurements including white spectral noise with sigma
    % std value.
    
    % extract parameters
    GM         = asterParams(1);
    Re         = asterParams(2);
    n_max      = asterParams(3);
    normalized = asterParams(4);
    
    W   = poleParams(1);
    W0  = poleParams(2);
    RA  = poleParams(3);
    DEC = poleParams(4);

    % output value 
    Nt  = length(time);
    ddU_ACI = ones(9, Nt);

    % time loop
    for j = 1:Nt
        rt_ACI = state(j, 1:3)';
        Wt = W0 + W * time(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        rt_ACAF = ACAF_ACI * rt_ACI;

        ACI_RTN = RTN2ECI(state(j, 1:3)', state(j, 4:6)');

        [~, ~, T_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rt_ACAF, Re, GM, normalized);
        T_ACI = ACAF_ACI' * T_ACAF * ACAF_ACI;
        
        ddU_ACI(:, j) = [T_ACI(1,1); T_ACI(1,2) ; T_ACI(1,3); T_ACI(2,1);...
         T_ACI(2,2); T_ACI(2,3) ; T_ACI(3,1); T_ACI(3,2); T_ACI(3,3)] + noise(:, j);

        % compute gravity partials
        [~, h_RTN] = potentialGradient_Cnm(n_max, rt_ACAF, Re, ...
            GM, (ACAF_ACI*ACI_RTN)', normalized);
        [~, h_ACI] = potentialGradient_Cnm(n_max, rt_ACAF, Re, ...
            GM, (ACAF_ACI)', normalized);
        H_ACI = h_ACI;
        H_RTN = h_RTN;
    end
end

