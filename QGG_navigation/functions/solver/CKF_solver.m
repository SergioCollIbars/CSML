function [X, Pt, Xhat, Xnot, pref, posf] = CKF_solver(TIME, X0, Xnot, P0, ...
    R0, Q0, Bw, T, planetParams, poleParams, C_mat, S_mat, ...
    system, consider_cov, augmented_st, DOM, posE, posM, posS)
    %%                    CKF FILTER FUNCTION
    % Description: Process measurements to refine real orbit. Use STM as  
    % a propagator.
    % Author: Sergio Coll Ibars
    % Date: 07/29/2024

    % compute dynamics. Use measured STM
    Nt = length(TIME);
    

    if(det(Bw) == 0)
        if(augmented_st), Ns = 7; else Ns = 6; end
        DMC = 0;  type = "SNC";
    else
        DMC = 1; Ns = 9; type = "DMC"; 
    end

    % state values
    X0 = X0 + Xnot;
    STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
    initState = [X0; STM0];

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, STATE] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, C_mat, S_mat, system, consider_cov, {DMC, Bw}, augmented_st), TIME, initState, options);
    X = STATE(:, 1:Ns)';
    STM  = STATE(:, Ns+1:end);

    Pt = ones(Nt, Ns*Ns) * NaN;
    Xhat = ones(Ns, Nt) * NaN;
    
    % data gap
    DG = 0;

    % start filter
    for j = 1:Nt
        % initial CKF and EKF estimations
        if(j == 1)
            X_hat = -Xnot;
            P = P0;
        end
        
        % compute previous STM
        if( j > 1)
            PHI_p = reshape(STM(j - 1, :), [Ns, Ns]);
            
            % STM @ ti 2 ti-1
            PHI_ij = reshape(STM(j, :), [Ns,Ns])/PHI_p;
        else
            % STM @ ti 2 ti-1
            PHI_ij = reshape(STM(j, :), [Ns,Ns]);
        end

         % compute gravity tensor measurements and partials
        state = X(:, j)';
        [Tc, Hx, ~] = compute_measurements(TIME(j), state, planetParams, ...
             poleParams, C_mat, S_mat, 0, consider_cov, augmented_st, [], DOM, posE(:, j), posM(:, j), posS(:, j), system);
        dY = T(:, j) - Tc;
        if(type == "DMC"), Hi = [Hx, zeros(6, 6)]; else,  Hi = Hx; end
         
        if(j == 1)
            Q = zeros(Ns,Ns);
            Qp = Q;
        else
            At = (TIME(j) - TIME(j-1)) / planetParams(3); % [sec]
            if(isnan(dY))
                DG = DG + At;
            else
                DG = 0;
            end
            At = At * planetParams(3);  % [-]
            DG = 0;
            Q = processNoise(Q0, DG, At, Bw, type, Ns);
        end
        Qp = Q;

        % run CKF
        [X_hat, P, ~, Qn] = CKF(dY, ...
             Hi, R0, P, X_hat, PHI_ij, Qp);
        Xhat(:, j) = X_hat;
    
        % addaptative process noise matrix
        if(type == "SNC" && det(Q0) ~= 0), Qp = Qn; end

        % update current state
        X(:, j) = X(:, j) + X_hat;

        % current uncertainty
        Pt(j, :) = reshape(P, [1, Ns*Ns]);
    end
    
    % restart filter
    PHI0 = reshape(STM(end, :), [Ns, Ns]);
    Xnot = PHI0\X_hat;

    % compute prefit and postfit
    Nm = length(R0);
    pref = ones(Nm, Nt) * NaN;
    posf = ones(Nm, Nt) * NaN;
    for j = 1:Nt
        PHI_i0 = reshape(STM(j, :), [Ns,Ns]);

        state = X(1:Ns, j)';
        [Tc, Hx, ~] = compute_measurements(TIME(j), state, planetParams,...
            poleParams, C_mat, S_mat, 0, 0, augmented_st, [], DOM, posE(:, j), posM(:, j), posS(:, j), system);
        if(type == "DMC"), Hi = [Hx, zeros(6, 6)]; else, ...
        Hi = Hx; end
        pref(:, j) = T(:, j) - Tc;
        posf(:, j) = (T(:, j) - Tc) - Hi * (PHI_i0 * Xnot);
    end
end

