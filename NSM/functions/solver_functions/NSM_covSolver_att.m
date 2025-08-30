function [sigma_N] = NSM_covSolver_att(planetParams, RotPlanet, B_ACI_mat, R, P0, Pc, Pxc, Xp, t ,...
    angVel, rn, vn, Iner, mask)
    % planet variables
     n_max = planetParams(3);

     % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);

     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % loop
    Ax_N = inv(P0); 
    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    fprintf('       Progress:    0%%');  % Initial message
    for j = 1:Nt-2
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);

        B_ACI = B_ACI_mat(minPos:maxPos, :);
   
        % Null space method correcting for attitude
        [Y_ACI, Hc_ACI, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp);
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];
        
        % compute attitude partials. Nominal body frame
        [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
        [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
        Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];

        % apply maks to partials
        hrot = Hrot(logical(mask), :);
        hc   =  Hc(logical(mask), :);

        C = null(hrot');
% %         [~,~,d] = svd(hrot');
% %         C = d(:, 3);

        hcp = C' * hc;
        r  = C' * R * C;
        hrotp = C' * hrot;
    
        ax = hcp' * (r\hcp);
        mxc = (hcp' * (r\hrotp));
        mcc = (hrotp' * (r\hrotp)); 
    
        Ax_N  = Ax_N + ax;
        Mxc = Mxc + mxc;
        Mcc = Mcc + mcc;

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

    disp(cond(Px));
    
    % compute correlations
    figure()
    R1 = P_N./ (sigma_N * sigma_N');
    figure;
    imagesc(abs(R1));
    colormap jet;
    colorbar;
    axis equal tight;
    title('Correlation Matrix. NSM');


    fprintf('\nDone!\n');
end

