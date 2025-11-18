function [X, Pc, Xhat, Xnot, pref, posf] = EKF_solver_EPHEM(TIME, X0, P0, ...
                    R0, Q0, Qb, meas, planetParams, BN_matrix, Cmat, Smat, ...
                    posE, posM, posS, gamma)
    % Run EKF solver using only the EPHEMERIDES dynamics.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % number of time steps
    Nt = length(TIME);
    Ns  = 12;

    % initiate filter
    STM0 = reshape(eye(Ns,Ns), [Ns*Ns,1]);
    X = zeros(Ns, Nt);
    X(:, 1) = X0;
    P = P0;

    % covariance & state correction at each time
    Pc = ones(Nt, Ns*Ns) * NaN;
    Xhat = ones(Ns, Nt) * NaN;
    
    % prefit and postfit
    Nm = length(R0);
    pref = ones(Nm, Nt) * NaN;
    posf = ones(Nm, Nt) * NaN;

    At = TIME(2) - TIME(1); % [-]

    delta_x = zeros(6, Nt);
    N = 10*Nt;

    % data gap null 
    f = waitbar(0, 'Starting');
    for k = 2:Nt
        % wait bar
        waitbar(k/Nt, f, sprintf('Progress: %d %%', floor(k/Nt*100)));

        % initial conditions
        t_span = [TIME(k - 1), TIME(k)];
        
        % Compute dynamics. ODE 113
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [~, STATE] = ode113(@(t, x) EOM_navigation_EPHEM(t, x,...
            planetParams, Cmat, Smat), t_span, [X0; STM0], options);
        
        state  = STATE(end, 1:Ns);
        STM    = STATE(:, Ns+1:end);
        PHI_ij = reshape(STM(end, :), [Ns,Ns]);

        % compute measurements
        maxInd = 3 * k; minInd = maxInd -2;
        BN = BN_matrix(minInd:maxInd, :);
        [T_c, ~] = compute_measurement_EPHEM(TIME(k), state, planetParams,...
            BN, Cmat, Smat, 0, posE(:, k), posM(:, k), posS(:, k));
        [Hmeas] = compute_meas_partials_EPHEM(TIME(k), state, planetParams,...
                    BN, Cmat, Smat, posE(:, k), posM(:, k), posS(:, k));

        % measurement residuals
        dY = meas(:, k) - T_c;
        pref(:, k) = dY;

% %         % distance to the Moon
% %         r_sc_moon = vecnorm(state(1:3)' - posM(:, k)).*planetParams(2);
% %         if(r_sc_moon < 1E7)
% %             Q_rv = Q0.*1E2; 
% %         else
% %             Q_rv = Q0;
% %         end
        Q_rv = Q0;

        % run filter
        Q = processNoise(Q_rv, 0, At, Qb, "SNC", Ns);
        
        if(k > N)
            Qxy = zeros(6, 6);
            for n = k-N:k-1
                Qxy = Qxy + delta_x(:, k-n) * delta_x(:, k-n)';
            end
            Qxy = Qxy./N;
            Qp = blkdiag(Qxy, Qb);
        else
            Qp = Q;
        end

        [X_hat, P, delta_x(:, k)] = EKF(dY, Hmeas, R0, P, PHI_ij, Qp);

        ee = sum(eig(P) < 0);
        if(ee > 0)
            disp(eig(P));
        end

        % update nominal
        X(:, k) = state' + X_hat;
        X0 = X(:, k);
        
        % postfit
        posf(:, k) = pref(:, k) - Hmeas * X_hat;

        % save covariance
        C = (P + P.')/2;
        P = C;
        Pc(k, :) = reshape(P, [1, Ns*Ns]);

        % save correction
        Xhat(:, k) = X_hat;
    end
    close(f);

    % end iterations
    Xnot = 1E-100*ones(Ns, 1);
end

