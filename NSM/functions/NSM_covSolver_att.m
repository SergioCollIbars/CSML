function [sigma_N] = NSM_covSolver_att(planetParams, RotPlanet, R, P0, Xp, t ,...
    attitude, datt_dt, ddatt_ddt, angVel, rn, vn)
    % pole & planet variables
    GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
    normalized = planetParams(4);  

     % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);
    
     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % loop
    Ax_N = inv(P0); 
    for j = 3:Nt-2
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
        [Hrot_dA_ang, Hrot_dAdT_ang] = compute_angularDyadPartials(angVel(:, j), attitude(:, j), datt_dt(:, j), ddatt_ddt(:, j));
        Hrot_dA = Hrot_grad + Hrot_dA_ang;
        Hrot_dAdT = Hrot_dAdT_ang;
        Hrot = [Hrot_dA, Hrot_dAdT];

        % Null space method correcting for attitude
        [~, ~, Hc] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp, B_ACI');
    
        C = null(Hrot');
        hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end); Hc(2, 2:end); Hc(5, 2:end);...
                Hc(8, 2:end);  Hc(3, 2:end); Hc(6, 2:end); Hc(9, 2:end)];
        hc = C' * hc;
        r  = C' * R * C;
    
        ax = hc' * inv(r) * hc;
    
        Ax_N  = Ax_N + ax;
    end

    Px = inv(Ax_N);
    sigma_N = sqrt(diag(Px));
end

