function [SH_N, sigma_N] = NSM_solver_att(planetParams, RotPlanet, R, P0, Xp, t ,...
    attitude, angVel, angAcc, Y, rn, vn, Iner)
    % pole & planet variables
    GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
    normalized = planetParams(4);  
    
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
% %             [Hrot_dA_ang, Hrot_dAdT_ang] = compute_angularDyadPartials(angVel(:, j), attitude(:, j), datt_dt(:, j), ddatt_ddt(:, j));
% %             Hrot_dA = Hrot_grad + Hrot_dA_ang;
% %             Hrot_dAdT = Hrot_dAdT_ang;
% %             Hrot = [Hrot_dA, Hrot_dAdT];

            [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
            Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];
        
            % Null space method correcting for attitude
            [Yc, ~, Hc_BODY] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                    zeros(9, Nt), Cp, Sp, B_ACI');
            [Yc] = add_angularComponents(Yc, attitude(:, j), zeros(3, Nt), angVel(:, j),...
                angAcc(:, j));
    
            [ax, nx] = NSM_method(Y(:,j), Yc, Hc_BODY, R, Hrot);

            Ax_N  = Ax_N + ax;
            Nx_N  = Nx_N + nx;
        end
    
        % solve LS
        XNOT_N = Ax_N\Nx_N;
    
        Xp(2:end) = Xp(2:end) + XNOT_N;
    
        [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);
    
        % update corrections
        xnot_N = xnot_N + XNOT_N;
    
        % show error
        disp('Null space update   '       + string(vecnorm(XNOT_N)));
        disp('--------------------------------------------------------'); 
        
        % update counter
        count = count + 1;
    end

    Px = inv(Ax_N);
    sigma_N = sqrt(diag(Px));
    
    [Xp_N] = mat2list(Cp, Sp, Nc, Ns);
    SH_N = Xp_N(2:end);
end


function [Hrot] = compute_rotPartials_analy(Y)
    Yxx = Y(1); Yxy = Y(2); Yxz = Y(3);
    Yyy = Y(5); Yyz = Y(6); Yzz = Y(9);

    Hrot = -[0, 2*Yxz, -2*Yxy;...
           -Yxz, Yyz, Yxx - Yyy;...
        Yxy, Yzz - Yxx, -Yyz;...
        -Yxz, Yyz, Yxx - Yyy;...
        -2*Yyz, 0, 2*Yxy;...
        Yyy - Yzz, -Yxy, Yxz;...
        Yxy, Yzz - Yxx, -Yyz;...
        Yyy - Yzz, -Yxy, Yxz;...
        2*Yyz, -2*Yxz, 0];
end

