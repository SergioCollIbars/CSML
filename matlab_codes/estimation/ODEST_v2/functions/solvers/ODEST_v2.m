function [Cnm_new, Snm_new, state, XNOT_c, XNOT_p, Pc_new, Pp_new, PHI_t] = ODEST_v2(t ,Y, rn, vn,...
        Cnm, Snm, poleParams, asterParams, sigma, XNOT_c, XNOT_p, ...
        P0_c, P0_p)
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

    % noise nominal & weights
    noise0 = zeros(9, Nt);
    scale = 5;
    Rc = diag([sigma, sigma, sigma, sigma, sigma])^2;
    Rp = diag([sigma * scale, sigma * scale, sigma * scale, sigma *scale].^2);
    
    % estimate gravity field
    Ax = inv(P0_c);
    Nx = -inv(P0_c) * XNOT_c;
    for j = 1:Nt
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % computed meas.
        [Yc, Hc, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cnm, Snm);

        % dY = [Y(2:3, j) - Yc(2:3);Y(5:6, j) - Yc(5:6)];
        dY = [Y(1:3, j) - Yc(1:3);Y(5:6, j) - Yc(5:6)];

        % sensitivity matrix positon partials & gravity partials
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, ...
            Re, GM, rn(:,j), ACAF_ACI);

        % look for null space
        C = null([Hpos(1:3, :); Hpos(5:6, :)]');

        % project measurements
        y  = C' * dY;
        hc = C' * [Hc(1:3, 2:end);Hc(5:6, 2:end)];
        r  = C' * Rc * C;

        % L.S filter
        Ax = Ax + (hc' * inv(r) * hc);
        Nx = Nx + (hc' * inv(r) * y);
    end
    XNOT_c = Ax\Nx;
    [X] = mat2list(Cnm, Snm, Nc, Ns);
    X(2:end) = X(2:end) + XNOT_c;
    Pc_new = inv(Ax);

    [Cnm_new, Snm_new] = list2mat(n_max, Nc, Ns, X);

    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [~, state] = ode113(@(t, x) EoM(t, x, Cnm_new, Snm_new, n_max, GM, Re, ...
    normalized, W0, W, RA, DEC), t, ...
    [rn(:, 1);vn(:, 1); reshape(eye(6,6), [36, 1])], options);

    rn = state(:, 1:3)';
    vn = state(:, 4:6)';
    PHI_t = state(:, 7:end);

    % solve for position
    Ax = inv(P0_p);
    Nx = -inv(P0_p) * XNOT_p;
    for j = 1:Nt
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

         % STM at current time
        PHI_i0 = reshape(PHI_t(j, :), [6,6]);
        
         % computed meas.
        [Yc, ~, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cnm_new, Snm_new);

        dY = [Y(4, j) - Yc(4); Y(7:9, j)- Yc(7:9)]; 
        
        [Hpos] = compute_posPartials(n_max, normalized, Cnm_new, Snm_new, ...
            Re, GM, rn(:,j), ACAF_ACI);
        Hvel = zeros(1, 3);
        Ht = [Hpos(4, :), Hvel;Hpos(7:9, :), zeros(3,3)];
        hp = Ht * PHI_i0;
        
        Ax = Ax + (hp' * inv(Rp) * hp);
        Nx = Nx + (hp' * inv(Rp) * dY);
    end
    XNOT_p = Ax\Nx;
    Pp_new = inv(Ax);

    state = ones(6, Nt);
    for j = 1:Nt
        PHI_i0 = reshape(PHI_t(j, :), [6,6]);
        state(:, j) = [rn(:, j);vn(:, j)] + PHI_i0 * XNOT_p;
    end
end

