function [X_smooth, P_smooth] = EKF_smoothing(time, planetParams, poleParams, Y, X0, P0, ...
    P0c_grav, BN_mat, Cnm, Snm, Q0, Qb, R0)
    
    Ns = length(X0);
    [~, ~, Ncs]  = count_num_coeff(planetParams(3)); 
    Nx = Ncs + Ns;

    % initiate filter
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    % output variables
    X_smooth = ones(Ns, length(time)) * NaN;
    P_smooth = nan(length(time), Ns*Ns);
    P_smooth(length(time), :) = reshape(P0, [Ns*Ns, 1]);

    Pxc = zeros(Ns, Ncs); Pcx = zeros(Ncs, Ns);
    Pcc = diag(diag(P0c_grav));
    P   = [P0, Pxc;...
           Pcx, Pcc];

    At = time(2) - time(1);

    for j = length(time):-1:2
            % initial conditions
            t_span = [time(j), time(j-1)];

            % current time index
            k = j - 1;
            
            % Compute dynamics. ODE 113
            [~, STATE] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
                Cnm, Snm), t_span, [X0;PHI0], options);
            
            state   = STATE(end, 1:Ns);
            PHI_tot = reshape(STATE(end, 13:end), [Nx, Nx]);
            
            % compute predicted measurements (TBD)
            maxInd      = 3 * k; minInd = maxInd - 2;
            BN          = BN_mat(minInd:maxInd, :);
            [Yc, Hi, Hc, ~]    = compute_measurements_filter(planetParams, ...
                poleParams, time(k), state, Cnm, Snm, BN);
            dY          = Y(:, k) - Yc; 
    
            % Include process noise
            Q = processNoise(Q0, Qb, At, Nx);
    
            % Run EKF process
            [X_hat, P] = EKF(dY, [Hi, Hc], R0, P, PHI_tot, Q, Ns, Ncs);
    
            ee = sum(eig(P) < 0);
            if(ee > 0)
                disp(eig(P));
            end
    
            % update nominal & save states
            X_smooth(:, k) = state' + X_hat(1:Ns);
            X0 = X_smooth(:, k);   

            P_smooth(k, :) = reshape(P(1:Ns, 1:Ns), [Ns*Ns, 1]);
    end
    
end

