function [X, P, Xhat, Xnot, pref, posf] = batch_solver(TIME, X0, Xnot, P0, ...
    R0, c, Pc, Pxc, T, planetParams, poleParams, count, C_mat, S_mat, ...
    system, consider_cov, augmented_st, DOM, posE, posM, posS)

    % compute dynamics. Use measured STM
    Nt   = length(TIME);
    Nc   = length(c);
    Ns = length(P0);
    Nm = length(T(:, 1));
    
    % state values
    X0 = X0 + Xnot;
    STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
    if(consider_cov) 
        SIGMA0 = reshape(zeros(6,Nc), [6*Nc, 1]); 
        initState = [X0; STM0; SIGMA0];
    else
        initState = [X0; STM0];
    end

    % integrate trajectory
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, STATE] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, C_mat, S_mat, system, consider_cov, {0,0}, augmented_st), TIME, ...
    initState, options);
    X = STATE(:, 1:Ns)';
    STM  = STATE(:, Ns+1:end);

    % consider parameters
    if(consider_cov)
        SIGMA = STATE(:, 43:end); 
        [Mxx_bar, Mxc_bar, Mcc_bar] = get_considerCov_apriori(P0, Pc, Pxc);
        Np = (Nc + 6)^2;
    else
        SIGMA =zeros(Nt, Ns*Nc); 
        Mxx_bar = inv(P0);
        Mxc_bar = zeros(Ns, Nc);
        Mcc_bar = zeros(Nc, Nc);
        c = c*0;
        Np = Ns*Ns;
    end
    
    P = ones(Nt, Np) * NaN;
    Xhat = ones(Ns, Nt) * NaN;
    Mxc = Mxc_bar;
    Mcc = Mcc_bar;
    
    % prefit and postfit
    pref = ones(Nm, Nt) * NaN;
    posf = ones(Nm, Nt) * NaN;

    % initialize information matrix & normal Eqn.
    if count == 0
        Ax = 0;
        Nx = 0;
    else
        Ax = Mxx_bar;
        Nx = - inv(P0) * Xnot;
    end

    % start filter
     [warnMsg, ~] = lastwarn;
     if isempty(warnMsg)
        for j = 1:Nt
            % compute previous STM and SIGMA
            PHI_i0 = reshape(STM(j, :), [Ns,Ns]);
            PSI_i0 = reshape(SIGMA(j, :), [Ns,Nc]);
            
            % compute correction
            Xhat(:, j) =  - PHI_i0 * Xnot;
    
            % compute gravity tensor measurements and partials
            state = X(:, j)';
            [Tc, Ht, ~] = compute_measurements(TIME(j), state, planetParams, ...
             poleParams, C_mat, S_mat, 0, consider_cov, augmented_st, [], DOM, posE(:, j), posM(:, j), posS(:, j), system);
           
            dY = T(:, j) - Tc;
            
            % relate visibility matrix to init time
            ht = Ht * PHI_i0;
            if(consider_cov), hc = Ht * PSI_i0 + Hc; else, hc = zeros(Nm, Nc); end
    
            % run batch
            if(isnan(dY))
                Ax = Ax + 0;
                Nx = Nx + 0;
            else
                Ax  = Ax  + (ht' * inv(R0) * ht);
                Nx  = Nx  + (ht' * inv(R0) * dY);
                Mxc = Mxc + (ht' * inv(R0) * hc);
                Mcc = Mcc + (hc' * inv(R0) * hc); 
            end
        end
        
        % compute initial deviation
        Xnot = Ax\Nx - Ax\(Mxc * c);
    
        % compute final covariance at epoch time
        Px = inv(Ax);
        Sxc = -Px * Mxc;
        Pxx = Px + Sxc*Pc*Sxc';
        Pxc = Sxc * Pc;
        
        if(consider_cov), P0 = [Pxx, Pxc;Pxc', Pc]; else, P0 = Px; end
        
        % compute prefit and postfit
        for j = 1:Nt
            PHI_i0 = reshape(STM(j, :), [Ns,Ns]);
            PSI_i0 = reshape(SIGMA(j, :), [Ns,Nc]);
            z = zeros(Nc,Ns);
            I = eye(Nc, Nc);
            if(consider_cov), PSI = [PHI_i0, PSI_i0; z, I]; else, PSI = PHI_i0; end
    
            p = PSI * P0 * PSI'; 
            P(j, :) = reshape(p, [1, Np]);
    
            state = X(:, j)';
            [Tc, Ht, ~] = compute_measurements(TIME(j), state, planetParams, ...
             poleParams, C_mat, S_mat, 0, consider_cov, augmented_st, [], DOM, posE(:, j), posM(:, j), posS(:, j), system);
            pref(:, j) = T(:, j) - Tc;
            posf(:, j) = (T(:, j) - Tc) - Ht * (PHI_i0 * Xnot);
        end
     end
end

