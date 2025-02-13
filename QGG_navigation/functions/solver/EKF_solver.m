function [X, Pc, Xhat, Xnot, pref, posf] = EKF_solver(TIME, X0, P0, ...
                    R0, Q0, Bw, T, planetParams, poleParams, Cmat, Smat, ...
                    system, consider_cov, augmented_st, DOM, posE, posM, posS)
    %%                    CKF FILTER FUNCTION
    % Description: Process measurements to refine real orbit.
    % Author: Sergio Coll Ibars
    % Date: 06/21/2024

    % number of time steps
    Nt = length(TIME);

    if(det(Bw) == 0)
        if(augmented_st), Ns = 7; else, Ns = 6; end
        DMC = 0;  type = "SNC";
    else
        DMC = 1; Ns = 9; type = "DMC"; 
    end

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

    % data gap null 
    DG = 0;
    f = waitbar(0, 'Starting');
    for k = 2:Nt
        % wait bar
        waitbar(k/Nt, f, sprintf('Progress: %d %%', floor(k/Nt*100)));

        % initial conditions
        t_span = [TIME(k - 1), TIME(k)];
        
        % Compute dynamics. ODE 113
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        [~, STATE] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
            poleParams, Cmat, Smat, system, 0, {DMC, Bw}, augmented_st), t_span, [X0; STM0], options);
        
        state = STATE(end, 1:Ns);
        STM  = STATE(:, Ns+1:end);
        PHI_ij = reshape(STM(end, :), [Ns,Ns]);
    
        [Tc, Hx, ~] = compute_measurements(TIME(k), state, planetParams, ...
             poleParams, Cmat, Smat, 0, consider_cov, augmented_st, [], DOM, posE(:, k), posM(:, k), posS(:, k), system);
        dY = T(:, k) - Tc; 
        pref(:, k) = dY;

        if(type == "DMC"), Hi = [Hx, zeros(6, 6)]; else, Hi = Hx; end

        % run filter
        At = (t_span(2) - t_span(1))/planetParams(3); % [sec]
% %         if(isnan(dY))
% %             DG = DG + At;
% %         else
% %             DG = 0;
% %         end
        At = At * planetParams(3);
        Q = processNoise(Q0, 0, At, Bw, type, Ns);
        Qp = Q; % don't use addaptative process noise

        if(k == 2), Qp = Q; end

        [X_hat, P, Qn] = EKF(dY, Hi, R0, P, PHI_ij, Qp);

        ee = sum(eig(P) < 0);
        if(ee > 0)
            disp(eig(P));
        end

        
        % addaptative process noise
        if(type == "SNC" && det(Q0) ~= 0), Qp = Qn; end

        % update nominal
        X(:, k) = state' + X_hat;
        X0 = X(:, k);
        
        % postfit
        posf(:, k) = pref(:, k) - Hi * X_hat;

        % save covariance
        Pc(k, :) = reshape(P, [1, Ns*Ns]);

        % save correction
        Xhat(:, k) = X_hat;
    end
    close(f);

    % end iterations
    Xnot = 1E-100*ones(Ns, 1);
end

