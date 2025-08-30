function [sigma_N, Sensitivity] = LS_covSolver_att(planetParams, RotPlanet, R, P0, Pc, Pxc, Xp, t ,...
    attitude, angVel, rn, vn, attErr, Iner, mask)
    % pole & planet variables
    GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
    normalized = planetParams(4); 

     % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);
    
     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    Ax_N = inv(P0);
    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    Serr = 0;
    fprintf('       Progress:    0%%');  % Initial message
    for j = 1:Nt-2
        % position vector
        rn_ACI = rn(:, j);
        
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);

        % from ACI to Nominal body frame
        B_ACI =rotationMatrix(attitude(1, j), attitude(2, j), attitude(3, j), ...
          [3, 2, 1]);   
        ACAF_B = ACAF_ACI * B_ACI';
    
         % compute attitude partials. Nominal body frame
        [Hrot_grad] = compute_rotPartials(n_max, normalized, Cp, Sp, Re, GM, rn_ACI, ACAF_ACI, ACAF_B);
        [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
        Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];
        
    
        % Least Square method
        [~, Hc_ACI] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp);
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);

        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];
        
        % apply maks to partials
        hrot = Hrot(logical(mask), :);
        hc   =  Hc(logical(mask), :);
    
        % information and normal matrices
        ax = hc' * inv(R) * hc;
        mxc = (hc' * inv(R) * hrot);
        mcc = (hrot' * inv(R) * hrot); 

        Ax_N  = Ax_N + ax;
        Mxc = Mxc + mxc;
        Mcc = Mcc + mcc;
        Serr = Serr + (hc' * inv(R) * hrot * attErr(1:6, j));

        % Update every ~5% (optional)
        if mod(j, round(Nt/20)) == 0 || j == Nt-2
            fprintf('\b\b\b\b%3d%%', round(100 * j / Nt));
        end
    end

    Px = inv(Ax_N);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pc*Sxc';
    Pxc = Sxc * Pc;
    P_N =  [Pxx, Pxc;Pxc', Pc];
    sigma_N = sqrt(diag(P_N));
    Sensitivity = - Px * Serr;
    
    disp(cond(Px));
    fprintf('\nDone!\n');
end

