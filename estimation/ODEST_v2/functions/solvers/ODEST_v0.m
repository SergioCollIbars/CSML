function [Cnm_new, Snm_new, state] = ODEST_v0(t ,Y, rn, vn,...
        Cnm, Snm, poleParams, asterParams, sigma, P0_c)
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
    Rc = diag([sigma, sigma])^2;
    
    % estimate gravity field
    xnot = zeros(Nc+Ns-1, 1);
    eps = 1E-10;
    err = 1;
    maxIter = 3;
    count = 0;
    while (err > eps) && (count < maxIter)
        Ax = inv(P0_c);
        Nx = -inv(P0_c) * xnot;
        for j = 1:Nt
            % RTN frame
            ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    
            % computed meas.
            [Yc, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                    noise0, Cnm, Snm);
    
            dy_RTN = ACI_RTN' * reshape(Y(:, j) - Yc, [3,3]) * ACI_RTN;
            dy = reshape(dy_RTN, [9, 1]);
            
            dY = [dy(5) - dy(9); dy(6)];
            hc = [Hc_RTN(5, 2:end) - Hc_RTN(9, 2:end); Hc_RTN(6, 2:end)];
    
            % L.S filter
            Ax = Ax + (hc' * inv(Rc) * hc);
            Nx = Nx + (hc' * inv(Rc) * dY);
        end
        XNOT_c = Ax\Nx;
        [X] = mat2list(Cnm, Snm, Nc, Ns);
         X(2:end) = X(2:end) + XNOT_c;
         [Cnm, Snm] = list2mat(n_max, Nc, Ns, X);

        % update states
        xnot = xnot + XNOT_c;
        err = vecnorm(XNOT_c);
        count = count + 1;
    end

    % estimate radial displacement from new Cnm / Snm coefficients
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    [~, state] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, ...
    normalized, W0, W, RA, DEC), t, ...
    [rn(:, 1);vn(:, 1); reshape(eye(6,6), [36, 1])], options);
    rn = state(:, 1:3)';
    vn = state(:, 4:6)';
    
    Ax = 0;
    Nx = 0;
    for j = 1:Nt
        % RTN frame
         ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    
        % computed meas.
        [Yc, ~, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cnm, Snm);

        dy_RTN = ACI_RTN' * reshape(Y(:, j) - Yc, [3,3]) * ACI_RTN;
        dy = reshape(dy_RTN, [9, 1]);
            
        dY = dy(1);
        H = -vecnorm(rn(:, j))^4 / (6*GM);

        Ax = Ax + (H' * inv(Rc(1,1)) * H);
        Nx = Nx + (H' * inv(Rc(1,1)) * dY);
    end
    
end


