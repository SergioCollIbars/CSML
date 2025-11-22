function [P0_new, state0_new] = initialize_filter(time, state0, Y, R0,...
    P0, Nt_max, planetParams, poleParams, Cnm, Snm, BN_mat)
    
    Ns = length(state0);
    PHI0 = reshape(eye(Ns,Ns), [Ns*Ns,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    Xnot = zeros(Ns, 1);
    err  = 1; thrs = 1E-15; maxIter = 20; count = 0;
    while((err > thrs) && (count < maxIter))
        % run Batch filter
        X0 = state0 + Xnot;
        [~, STATE] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
            Cnm, Snm), time(1:Nt_max), [X0;PHI0], options);

        X    = STATE(:, 1:Ns)';
        STM  = STATE(:, Ns+1:end);

        Ax = inv(P0); Nx = -P0\Xnot;
        for k = 1:Nt_max
            PHI_i0 = reshape(STM(k, :), [Ns,Ns]);

            % compute measurements and partials
            maxInd      = 3 * k; minInd = maxInd - 2;
            BN          = BN_mat(minInd:maxInd, :);
            [Yc, Hi]    = compute_measurements_filter(planetParams, ...
                poleParams, time(k), X(:, k)', Cnm, Snm, BN);
            dY          = Y(:, k) - Yc;
            H0          = Hi * PHI_i0;

            Ax  = Ax  + (H0' * (R0\ H0));
            Nx  = Nx  + (H0' * (R0\dY));
        end
        
        Xhat   = Ax\Nx;
        Xnot   = Xnot + Xhat;

        err    = vecnorm(Xhat); disp('  err: ' + string(err));
        count  = count + 1;
    end

    P0_new     = inv(Ax);
    state0_new = X(:, 1);
end

