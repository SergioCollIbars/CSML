function [SH_N, sigma_N] = LS_solver_att(planetParams, RotPlanet, B_ACI_mat, R, P0, Pc, Pxc, Xp, t ,...
    attitude, angVel, angAcc, Y, rn, vn, Iner, mask)
    % planet variables
    n_max = planetParams(3);
 
    
    % variables number
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
     Nt          = length(t);
    
     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % loop
    iterMax = 3;
    count   = 0;
    xnot_N = zeros(Ncs-1, 1);
    while count < iterMax
        Ax_N = inv(P0); Nx_N = -inv(P0) * xnot_N;
        [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
        for j = 1:Nt-2
            % Planet orientation
            maxPos = 3*j; minPos = maxPos - 2;
            ACAF_ACI = RotPlanet(minPos:maxPos, :);
    
            % from ACI to Nominal body frame
            B_ACI = B_ACI_mat(minPos:maxPos, :);
% %             ACAF_B = ACAF_ACI * B_ACI';
        
             % compute attitude partials. Nominal body frame
% %             [Hrot_grad] = compute_rotPartials(n_max, normalized, Cp, Sp, Re, GM, rn_ACI, ACAF_ACI, ACAF_B);
% %             [Hrot_dA_ang, Hrot_dAdT_ang] = compute_angularDyadPartials(angVel(:, j), attitude(:, j), datt_dt(:, j), ddatt_ddt(:, j), Iner);
% %             Hrot_dA = Hrot_grad + Hrot_dA_ang;
% %             Hrot_dAdT = Hrot_dAdT_ang;
% %             Hrot = [Hrot_dA, Hrot_dAdT];
        
            % Null space method correcting for attitude
            [Y_ACI, Hc_ACI, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                    zeros(9, Nt), Cp, Sp);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
                  Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];

            [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
            [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
            Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];

            [Yc] = add_angularComponents(Y_ACI, attitude(:, j), zeros(3, Nt), angVel(:, j),...
                angAcc(:, j));
            
            [ax, nx, mxc, mcc] = LS_method(Y(:,j), Yc(logical(mask)),...
                Hc(logical(mask), :), Hrot(logical(mask), :), R);

            Ax_N  = Ax_N + ax;
            Nx_N  = Nx_N + nx;
            Mxc = Mxc + mxc;
            Mcc = Mcc + mcc;
        end
    
        % solve LS
        XNOT_N = Ax_N\Nx_N;
    
        Xp(2:end) = Xp(2:end) + XNOT_N;
    
        [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);
    
        % update corrections
        xnot_N = xnot_N + XNOT_N;
    
        % show error
        disp('  Least Squares update   '       + string(vecnorm(XNOT_N)));
        disp('  Condt. Numb. =    '          + string(cond(Ax_N)));
        disp('--------------------------------------------------------'); 
        
        % update counter
        count = count + 1;
    end

    Px = compute_covarianceMat(Ax_N);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pc*Sxc';
    Pxc = Sxc * Pc;
    P_N =  [Pxx, Pxc;Pxc', Pc];
    
    sigma_N = sqrt(diag(P_N));
    
    [Xp_N] = mat2list(Cp, Sp, Nc, Ns);
    SH_N = Xp_N(2:end);
end
