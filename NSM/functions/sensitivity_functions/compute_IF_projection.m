function [IF_proj,IF_init, cond_proj, cond_init] = compute_IF_projection(planetParams, RotPlanet, ...
    B_ACI_mat, R, invP0, Xp, t, angVel, rn, vn, Iner, mask)
    % COMPUTE INFORMATION MATRIX PROJECTION
    % Description: Given a trajectory and orientation data set, compute the
    % projection of the IF in the projected and original space.
    % Author: Sergio Coll Ibars
    % Date: 06/29/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % planet variables
    n_max = planetParams(3);
    
    % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);
    
     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % store condition number over time
    cond_proj = zeros(1, Nt); cond_init = zeros(1, Nt);
    Ax_N = invP0; Ax = invP0;
    fprintf('       Progress:    0%%');  % Initial message
    for j = 3:Nt-2
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);

        % from ACI to Nominal body frame
        B_ACI = B_ACI_mat(minPos:maxPos, :);

        % angular vel / acc dyads partials
        [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
       
        % Null space method correcting for attitude
        [Y_ACI, Hc_ACI, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp);

        % rotate coefficient partials to body frame
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
              Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];
        
        %  tensor rotation partials
        [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);

        % ensamble partials
        Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];

        % apply maks to partials
        hrot = Hrot(logical(mask), :);
        hc   =  Hc(logical(mask), :);

        % apply projection
        C = null(hrot');
        r = C' * R * C;

        % projected IF & condition number
        Ax_N  = Ax_N + (C' * hc)' * (r\(C' * hc));
        % cond_proj(j) = cond(Ax_N);

        % original IF
        Ax = Ax  +  hc' * (R\hc);
        % cond_init(j) = cond(Ax);

        % Update every ~5% (optional)
        if mod(j, round(Nt/20)) == 0 || j == Nt-2
            fprintf('\b\b\b\b%3d%%', round(100 * j / Nt));
        end
    end
    IF_proj = Ax_N; IF_init = Ax;
    fprintf('\nDone!\n');
end

