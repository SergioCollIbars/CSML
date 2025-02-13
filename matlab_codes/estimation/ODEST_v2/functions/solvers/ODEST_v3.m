function [Cnm_new, Snm_new, state, Pc_new, Pp_new, PHI_t] = ODEST_v3(t ,Y, rn, vn,...
        Cnm, Snm, poleParams, asterParams, sigma, P0_c, P0_p)
    %ODEST_V2 estimation process based on null space
    %   Estimate position and gravity field based on null space of sensitivity
    %   matrix
    
    % extract parameters
    GM         = asterParams(1);
    Re         = asterParams(2);
    n_max      = asterParams(3);
    normalized = asterParams(4);
    [Nc, Ns, ~] = count_num_coeff(n_max);

    W   = poleParams(1);
    W0  = poleParams(2);
    RA  = poleParams(3);
    DEC = poleParams(4);
    
    % output value 
    Nt = length(t);

    % noise nominal
    noise0 = zeros(6, Nt);

    % max iterations
    MaxIter = 10;
    
    % estimate gravity field
    xnot_c = zeros(Nc + Ns - 1, 1);
    eps = 1E-15;
    err = eps + 1;
    count = 0;
    while (err > eps) && (count < MaxIter)
        Ax = inv(P0_c);
        Nx = -inv(P0_c) * xnot_c;
        Rc = diag(sigma)^2;
        for j = 1:Nt
            % computed meas.
            [Yc, Hc] = gradiometer_meas(t(j) ,asterParams, poleParams, rn(:, j)', ...
                    noise0, Cnm, Snm);
    
            dY = Y(1:4, j) - Yc(1:4);
    
            % sensitivity matrix positon partials & gravity partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, ...
                Re, GM, rn(:,j));
    
            % look for null space
            C = null(Hpos(1:4, :)');
    
            % project measurements
            y = C' * dY;
            hc = C' * Hc(1:4, 2:end);
            r = Rc;
    
            % L.S filter
            Ax = Ax + (hc' * inv(r) * hc);
            Nx = Nx + (hc' * inv(r) * y);
        end
        XNOT_c = Ax\Nx;
        [X] = mat2list(Cnm, Snm, Nc, Ns);
        X(2:end) = X(2:end) + XNOT_c;

        % update gravity field
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X);

        % update correction
        xnot_c = xnot_c + XNOT_c;
        err = vecnorm(XNOT_c);
        disp('error SH: ' + string(err) + 'Iteration = ' + string(count))

        % update counter 
        count = count + 1;
    end
    
    % new grav field to estimate position errors
    [Cnm_new, Snm_new] = list2mat(n_max, Nc, Ns, X);
    Pc_new = inv(Ax);

    % solve for position
    xnot_p = zeros(6, 1);
    eps = 1E-10;
    err = eps + 1;
    count = 0;
    while (err > eps) && (count < MaxIter)
        % new trajectory 
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        [~, state] = ode113(@(t, x) EoM(t, x, Cnm_new, Snm_new, n_max, GM, Re, ...
        normalized, W0, W, RA, DEC), t, ...
        [rn(:, 1);vn(:, 1); reshape(eye(6,6), [36, 1])], options);
    
        rn = state(:, 1:3)';
        vn = state(:, 4:6)';
        PHI_t = state(:, 7:end);

        Ax = inv(P0_p);
        Nx = -inv(P0_p) * xnot_p;
        Rp = diag([sigma, sigma].^2);
        for j = 1:Nt
             % STM at current time
            PHI_i0 = reshape(PHI_t(j, :), [6,6]);
            
             % computed meas.
            [Yc, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, rn(:, j)', ...
                    noise0, Cnm_new, Snm_new);
            
            dY = Y(5:6, j) - Yc(5:6);
            
            [Hpos] = compute_posPartials(n_max, normalized, Cnm_new, Snm_new, ...
                Re, GM, rn(:,j));
            Hvel = zeros(2, 3);
            Ht = [Hpos(5:6, :), Hvel];
            hp = Ht * PHI_i0;
            
            Ax = Ax + (hp' * inv(Rp) * hp);
            Nx = Nx + (hp' * inv(Rp) * dY);
        end
        XNOT_p = Ax\Nx;
        
        % update initial condition
        rn(:, 1) = rn(:, 1) + XNOT_p(1:3);
        vn(:, 1) = vn(:, 1) + XNOT_p(4:6);

        % update correction
        xnot_p = xnot_p + XNOT_p;
        err = vecnorm(XNOT_p(1:3));
        disp('error Pos: ' + string(err) + 'Iteration = ' + string(count))

        % update counter 
        count = count + 1;
    end

    Pp_new = inv(Ax);

    state = ones(6, Nt) * NaN;
    for j = 1:Nt
        PHI_i0 = reshape(PHI_t(j, :), [6,6]);
        state(:, j) = [rn(:, j);vn(:, j)] + PHI_i0 * XNOT_p;
    end
end

