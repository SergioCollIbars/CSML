function [X, Pc, Xhat, Xnot, pref, posf] = UKF_solver(time, X0, P0, ...
                    R0, Q0, Bw, T, planetParams, poleParams, Cmat, Smat, ...
                    system, consider_cov, augmented_st,DOM, posE, posM, posS)
    %%                    CKF FILTER FUNCTION
    % Description: Process measurements to refine real orbit.
    % Author: Sergio Coll Ibars
    % Date: 06/21/2024
    
    % number of time steps
    Nt = length(time);
    Nm = length(T(:, 1));

    if(det(Bw) == 0), DMC = 0; Ns = 6; type = "SNC"; ...
    else, DMC = 1; Ns = 9; type = "DMC"; end

    % initiate filter
    X = zeros(Ns, Nt);
    X(:, 1) = X0;
    STM0 = reshape(eye(Ns, Ns), [Ns*Ns, 1]);

    % covariance & state correction at each time
    Pc = ones(Nt, Ns*Ns) * NaN; Pc(1, :) = reshape(P0, [Ns*Ns, 1]);
    Xhat = ones(Ns, Nt) * NaN;
    
    % UKF routine
    DG = 0;
    f = waitbar(0, 'Starting');
    for k = 2:Nt
        % wait bar
        waitbar(k/Nt, f, sprintf('Progress: %d %%', floor(k/Nt*100)));

        % create sigma points. State @ k-1
        Xhat_prev = X(:, k-1);
        [xhat_i_prev] = sigmaPoints_state(Ns, Pc(k-1, :), Xhat_prev);
    
        % propagate sigma points using N.L funcitons. State @ k
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        t_span = [time(k-1), time(k)];
        Xhat_i = xhat_i_prev.*0;
        for i = 1:2*Ns   
            [~, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
            poleParams, Cmat, Smat, system, 0, {DMC, Bw}, augmented_st), t_span, [xhat_i_prev(:, i); STM0], options);
            Xhat_i(:, i) = state(end, 1:6)';
        end
    
        % apriori state estimate @ time k
        Xhat_min = 1/(2*Ns).* sum(Xhat_i, 2);
    
        % apriori state covariance @ time k
        A  = zeros(Ns, Ns);
        for i = 1:2*Ns
            A = A + (Xhat_i(:, i) - Xhat_min) * (Xhat_i(:, i) - Xhat_min)';
        end
        At = time(k) - time(k-1);
        Q = processNoise(Q0, DG, At, Bw, type, Ns);
        P_min = 1/(2*Ns).* A + Q;
    
         % create sigma points. State @ k
        [Xhat_i] = sigmaPoints_state(Ns, reshape(P_min, [36, 1]), Xhat_min);
    
        % compute predicted measurements. State @ time k
        Yhat_i = zeros(Nm, 2*Ns);
        for i = 1:2*Ns
            [Yhat_i(:, i), ~, ~] = compute_measurements(time(k), Xhat_i(1:6,i)', planetParams, ...
             poleParams, Cmat, Smat, 0, consider_cov, augmented_st, [], DOM, posE(:, k), posM(:, k), posS(:, k), system);
        end
        Yhat = 1/(2*Ns).* sum(Yhat_i, 2);
        
        if isnan(vecnorm(T(:, k)))
            X(:, k) = Xhat_min;
            Pc(k, :) = reshape(P_min, [36, 1]);
        else
            % compute measurement covariance
             B  = zeros(Nm, Nm);
            for i = 1:2*Ns
                B = B + (Yhat_i(:, i) - Yhat) * (Yhat_i(:, i) - Yhat)';
            end
             Py = 1/(2*Ns).* B + R0;
             
            C  = zeros(Ns, Nm);
            for i = 1:2*Ns
                C = C + (Xhat_i(:, i) - Xhat_min) * (Yhat_i(:, i) - Yhat)';
            end
             Pxy = 1/(2*Ns).* C;
        
             % kalman update
             K = Pxy/(Py);
             Xhat_plus = Xhat_min + K * (T(:, k) - Yhat);
             P_plus = P_min - K * Py * K';

            ee = sum(eig(P_plus) < 0);
            if(ee > 0)
                s = 1;
            end
            
             % save states
             X(:, k) = Xhat_plus;
             Pc(k, :) = reshape(P_plus, [36, 1]);
        end
    end
    close(f);

    % end iterations
    Xnot = 1E-100*ones(Ns, 1);

    pref = ones(6, Nt) * NaN;
    posf = ones(6, Nt) * NaN;
end



function [xhat_i] = sigmaPoints_state(Ns, P, Xhat)
    xhat_i = zeros(Ns, 2*Ns);
    xhat_i(:, :) = Xhat.* ones(Ns, 2*Ns);
    % compute matrix square root
    L = chol(Ns.*reshape(P, [6,6]), 'lower');
    X = L';

    % propagate first set of points
    for i = 1:Ns
        xhat_i(:, i)    =  xhat_i(:, i)    + X(i, :)';
        xhat_i(:, i+Ns) =  xhat_i(:, i+Ns) - X(i, :)';
    end
end

