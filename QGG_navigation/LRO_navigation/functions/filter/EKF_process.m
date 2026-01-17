function [X_EKF, P_EKF, posfit] = EKF_process(time, planetParams, Y_N, ...
    signal_error, noise_att, offset_att, sigma_att, X0_N, P0_N, ...
    orientation, NB_EARTH, NB_MOON, Cnm_list, Snm_list, Q0, Qb, R0, mask)
    
    % state mask (pos, vel & bias)
    mask_state = [ones(6, 1);mask];

    % Number of total states
    Nx = length(X0_N);
    
    % Number of trajectory states (pos & vel)
    Ns = 6;

    % initiate filter
    PHI0 = reshape(eye(Ns,Ns), [Ns*Ns,1]);
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);

    if(orientation == "Inertial")
        BN0 = eye(3);
    elseif(orientation == "RTN")            
        NB = RTN2ECI(X0_N(1:3), X0_N(4:6));
        BN0= NB';
    end

    % output variables
    X_EKF = ones(Nx, length(time)) * NaN;
    P_EKF = nan(length(time), Nx*Nx);
    
    % rotate initial uncertainty to body frame
    P0          = rotate_P0(BN0,P0_N);
    P_old       = P0;
    Pc          = eye(3).*(sigma_att)^2;

    % assign initial values
    P_EKF(1, :) = reshape(P0_N, [Nx*Nx, 1]);
    X_EKF(:, 1) =  X0_N;

    % postfit measurements
    posfit  = nan(6, length(time));
    
    At = time(2) - time(1);
    fprintf('  0%%');
    for k = 2:length(time)
        pct = 100 * k / length(time);
        fprintf('\b\b\b\b%3.0f%%', pct);
    
        % initial conditions
        t_span = [time(k-1), time(k)];
        
        % Compute dynamics. ODE 113 (Inertial frame)
        [~, STATE] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
            Cnm_list, Snm_list), t_span, [X0_N(1:Ns);PHI0], options);
        
        state_N   = STATE(end, 1:Ns);
        PHI_N     = reshape(STATE(end, Ns+1:end), [Ns, Ns]);
        
        % compute S/C orientation 
        [BN_mat] = compute_orientation_SC(t_span, ...
            [STATE(1, 1:Ns);STATE(end, 1:Ns)], orientation);

        maxInd         = 3 * k; minInd = maxInd - 2;
        BN0            = BN_mat(1:3, :);
        BN             = BN_mat(4:6, :);
        BODYMOON_J2000 = NB_MOON(minInd:maxInd, :)';

        % cumpute true measurements in S/C body frame
        BT2_BT1     = rotationMatrix(noise_att(1, k), ...
                    noise_att(2, k), noise_att(3, k), [1, 2, 3]); 
        BT1_B       = rotationMatrix(offset_att(1, k), ...
                    offset_att(2, k), offset_att(3, k), [1, 2, 3]); 
        [Y]      = compute_orientation_Meas(time(k), BT2_BT1*BT1_B*BN, Y_N(:, k), ...
                   signal_error(:, k));
        
        % rotate to body frame
        state     = rotate_state(time(k), BN, state_N');
        PHI_vec   = rotate_STM(time(k), BN, reshape(PHI_N, [1, Ns*Ns]), BN0);
        PHI       = reshape(PHI_vec, [Ns, Ns]);

        PHI_tot   = [PHI, zeros(6,6);zeros(6,6), eye(6)];
        
        % actual RTN frame
        [Yc, Hi, Hrot_i]  = compute_measurements_filter(planetParams, time(k), ...
                       state_N, Cnm_list, Snm_list, BN, BODYMOON_J2000);
        dY           = Y - Yc - X0_N(Ns+1:end); 
    
        % Include process noise
        Q = processNoise(Q0, Qb, At, Nx);
        % % Q  = rotate_processNoise(BN,Q_N, 1);

        % apply measurement mask
        dY_used = dY(logical(mask));
        Hi_used = Hi(logical(mask), logical(mask_state));
        Hrot_i_used = Hrot_i(logical(mask), :);
        R0_used = R0(logical(mask), logical(mask));

        % consider covarinace inflation
        R_inflation = Hrot_i_used * Pc * Hrot_i_used';
        R0_used     = R0_used + R_inflation;

        % apply state mask
        PHI_used = PHI_tot(logical(mask_state), logical(mask_state));
        Q_used   = Q(logical(mask_state), logical(mask_state));
        P_used   = P_old(logical(mask_state), logical(mask_state));
    
        % Run EKF process
        [X_hat_new, P_new] = EKF(dY_used,Hi_used, R0_used, P_used, ...
            PHI_used, Q_used);
    
% %         ee = sum(eig(P_new) < 0);
% %         if(ee > 0)
% %             disp(eig(P_new));
% %         end
        
        % update states inside mask (body frame)
        X_hat = zeros(Nx, 1); X_hat(logical(mask_state)) = X_hat_new;
        P_old(logical(mask_state), logical(mask_state))  = P_new;
                
        % update nominal & save states (Inertial frame)
        A       = blkdiag(BN', BN', eye(6));
        XB      = state + X_hat(1:Ns);
        XB_bias = X0_N(Ns+1:end) + X_hat(Ns+1:end);

        X_EKF(:, k)     = A * [XB;XB_bias];

        % compute posfit
        posfit(:, k) = dY - Hi * X_hat;

        % rotate to inertial state
        X0_N = X_EKF(:, k);   
        
        P_N  = A * P_old * A.';
        P_EKF(k, :) = reshape(P_N, [Nx*Nx, 1]);

        % update Body frame definition
        [BN_mat] = compute_orientation_SC(t_span, ...
            X_EKF(1:Ns, k-1:k)', orientation);
        BN     = BN_mat(4:6, :); 
        A      = blkdiag(BN, BN, eye(6));
        P_old  = A * P_N * A.'; 
    end
    fprintf('\n');
end

