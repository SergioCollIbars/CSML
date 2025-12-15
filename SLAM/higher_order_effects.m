function [ddU_ACI] = higher_order_effects(time ,asterParams, RotPlanet, state, ...
                noise, Cnm, Snm, err)

    % extract parameters
    GM         = asterParams(1);
    Re         = asterParams(2);
    n_max      = asterParams(3);
    normalized = asterParams(4);

    % output value 
    Nt  = length(time);
    ddU_ACI = ones(9, Nt);

    if(Nt == 1)
        j = 1;
        rt_ACI = state(j, 1:3)';
        
        % rotation matrix
        maxPos = 3* j; minPos = maxPos -2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);
        rt_ACAF = ACAF_ACI * rt_ACI;

        [~, ~, T_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rt_ACAF, Re, GM, normalized);
        T_ACI = ACAF_ACI' * T_ACAF * ACAF_ACI;
        
        ddU_ACI(:, j) = [T_ACI(1,1); T_ACI(1,2); T_ACI(1,3); T_ACI(2,1);...
         T_ACI(2,2); T_ACI(2,3); T_ACI(3,1); T_ACI(3,2); T_ACI(3,3)] + noise(:, j);
        
        % compute gravity partials
        [~, h_ACI] = potentialGradient_Cnm(n_max, rt_ACAF, Re, ...
            GM, (ACAF_ACI)', normalized);
        H_ACI = h_ACI;
        H_BODY = [];
    else
        % time loop
        for j = 1:Nt
            rt_ACI = state(j, 1:3)';
            
            % rotation matrix
            maxPos = 3* j; minPos = maxPos -2;
            ACAF_ACI = RotPlanet(minPos:maxPos, :);
            rt_ACAF = ACAF_ACI * rt_ACI;
    
            [~, ~, T_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                               rt_ACAF, Re, GM, normalized);
            T_ACI = ACAF_ACI' * T_ACAF * ACAF_ACI;
            
            ddU_ACI(:, j) = [T_ACI(1,1); T_ACI(1,2) ; T_ACI(1,3); T_ACI(2,1);...
             T_ACI(2,2); T_ACI(2,3) ; T_ACI(3,1); T_ACI(3,2); T_ACI(3,3)] + noise(:, j);

            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                    rt_ACI, ACAF_ACI, ACAF_ACI);
            
            first_error = Hpos * err(:, j);

            ddU_ACI(:, j) = ddU_ACI(:, j) - first_error;
        end
        % empty arrays
        H_ACI = [];
        H_BODY = [];
    end

end

