function [] = compute_eigenvalue(planetParams, RotPlanet, B_ACI_mat, ...
        Xp, t, angVel, rn, vn, Iner)
     % planet variables
     n_max = planetParams(3);

     % variables number
    [Nc, Ns, ~] = count_num_coeff(n_max); 
     Nt          = length(t);

     % a-priori coefficients
    [Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

    % output values
    eig_full  = zeros(9, Nt);
    eig_six   = zeros(6, Nt); 
    eig_three = zeros(3, Nt);

    % Time normalization (meas. sampling time)
    T = 10;                 % [s]

    fprintf('       Progress:    0%%');  % Initial message
    for j = 1:Nt-2
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);

        B_ACI = B_ACI_mat(minPos:maxPos, :);
   
        % Null space method correcting for attitude
        [Y_ACI, ~, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                zeros(9, Nt), Cp, Sp);

        % compute attitude partials. Nominal body frame
        [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(angVel(:, j), Iner);
        [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
        hrot = [Hrot_grad, (Hrot_omega_dyad+H_omegaDot_dyad)./T]./(1E-9);   % [E/rad]
        
        % compute Eigenvalue full tensor
        [~, v, ~] = svd(hrot');
        eig_full(:, j) = sum(v)';

        % compute Eigenvalue 6 components
        mask       =  [1, 1, 1, 0, 1, 1, 0, 0, 1]';
        [~, v, ~] = svd(hrot(logical(mask), :)');
        eig_six(:, j) = sum(v)';

        % compute Eigenvalue 3 components
        mask       =  [1, 0, 0, 0, 1, 0, 0, 0, 1]';
        [~, v, ~] = svd(hrot(logical(mask), :)');
        eig_three(:, j) = sum(v)';

         % Update every ~5% (optional)
        if mod(j, round(Nt/20)) == 0 || j == Nt-2
            fprintf('\b\b\b\b%3d%%', round(100 * j / Nt));
        end
    end
    
    gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
    t_UTC = gps_epoch + seconds(t);        % date time 
    t = t_UTC;
         
    figure()
    semilogy(t, eig_full + 1E-20, 'LineWidth', 2)
    title('Eigenvalues full tensor')
    grid on;
    ylabel('[Eotvos / rad]')

    figure()
    semilogy(t, eig_six, 'LineWidth', 2)
    title('Eigenvalues six components')
    grid on;
    ylabel('[Eotvos / rad]')

    figure()
    semilogy(t, eig_three, 'LineWidth', 2)
    title('Eigenvalues three components')
    grid on;
    ylabel('[Eotvos / rad]')

end

