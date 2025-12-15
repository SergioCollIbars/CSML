function [P0_new, state0_new] = initialize_filter(time, state0, Y, R0,...
    P0, Pc, Nt_max, planetParams, poleParams, Cnm, Snm, BN_mat)
    
    Ns = length(state0);
    [~, ~, Ncs]  = count_num_coeff(planetParams(3)); 
    Nxx = Ncs + Ns;

    PHI0 = reshape(eye(Nxx,Nxx), [Nxx*Nxx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % consider parameters
    % % Pcc = diag(diag(Pc));
    Pcc = Pc;

    Xnot = zeros(Ns, 1);
    err  = 1; thrs = 1E-15; maxIter = 10; count = 0;
    Mx = zeros(Ns, Ncs);
    while((err > thrs) && (count < maxIter))
        % run Batch filter
        X0 = state0 + Xnot;
        [~, STATE] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
            Cnm, Snm), time(1:Nt_max), [X0;PHI0], options);

        X    = STATE(:, 1:Ns)';
        STM  = STATE(:, Ns+1:end);

        Ax = inv(P0); Nx = -P0\Xnot;
        for k = 1:Nt_max
            PHI_tot = reshape(STM(k, :), [Nxx,Nxx]);
            PHI_x   = PHI_tot(1:Ns, 1:Ns);
            PHI_xc  = PHI_tot(1:Ns, Ns+1:end);

            % compute measurements and partials
            maxInd      = 3 * k; minInd = maxInd - 2;
            BN          = BN_mat(minInd:maxInd, :);
            [Yc, Hi, Hc, ~]    = compute_measurements_filter(planetParams, ...
                poleParams, time(k), X(:, k)', Cnm, Snm, BN);
            dY          = Y(:, k) - Yc;
            H0          = Hi * PHI_x;
            H0c         = Hc + Hi * PHI_xc;

            Ax  = Ax  + (H0' * (R0\ H0));
            Nx  = Nx  + (H0' * (R0\dY));

            Mx  = Mx + H0' * (R0\H0c);
        end
        
        Xhat   = Ax\Nx;
        Xnot   = Xnot + Xhat;

        err    = vecnorm(Xhat);
        count  = count + 1;
    end
    P_xx        = inv(Ax);
    Sxc         = -P_xx * Mx;
    Pxc         = Sxc * Pcc;
    P0_new      = [P_xx+Sxc*Pcc*Sxc',Pxc;Pxc',Pcc];

    state0_new = X0;
end

