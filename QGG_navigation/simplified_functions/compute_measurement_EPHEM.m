function [T, acc] = compute_measurement_EPHEM(t, state, planetParams,...
    BN_matrix, C_mat, S_mat, Nn, posE, posM, posS)
    % Compute measurements using the EPHEMERIDES model. The measurements
    % computed are: gradiometer measurements and acceloremeter measurements
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % S/C paramters
    mass = planetParams(11);
    A    = planetParams(12);

    % Simulation points
    Nt  = length(t);
    
    % output matrices
    T   = ones(6, Nt) * NaN;
    acc = ones(3, Nt) * NaN;

    % Generate noise
    measDim_QGG = planetParams(3)^2;
    measDim_Acc = planetParams(2)*planetParams(3)^2;
    Hz = planetParams(3);
    if(Nt > 1)
        dt = (t(2) - t(1)) / planetParams(3);                        % [sec]
        frec = 1 / dt;
        At = t(2) - t(1);                                            % [-]
    else
        frec = 1;
        At = 1;
    end
    sigma_GG  = [1, 1, 1, 1, 1, 1] * 1E-12 * sqrt(frec);             % [1/s^2]
    sigma_A   = 1E-10 * sqrt(frec);                                  % [m/s^2]
    sigma_bGG = ones(1, 6) * 9.33E-18;                               % [1/s^2 * sqrt(Hz)] works with 9E-18
    sigma_bac = ones(1, 3) * 9.33E-16;                                      % [m/s^2 * sqrt(Hz)]
    sigma_GG  = sigma_GG./measDim_QGG;                               % [-]
    sigma_A   = sigma_A./measDim_Acc;                                % [-]
    sigma_bGG = sigma_bGG./(measDim_QGG * sqrt(Hz));                 % [-]
    sigma_bac = sigma_bac./(measDim_Acc * sqrt(Hz));                 % [-]

    
    noise_GG = ones(6, Nt);
    for j=1:length(sigma_GG)
            noise_GG(j, :) = normrnd(0, sigma_GG(j), [1, Nt]) * Nn;
    end
    noise_A = normrnd(0, sigma_A, [3, Nt]) * Nn;

    % include bias random-walk
    if(Nn == 1)
        for j = 1:Nt-1
            for k = 1:length(sigma_bGG)
                state(j+1, 7+k) = state(j, 7+k) + normrnd(0,...
                    sqrt(At)*sigma_bGG(k), [1, 1]);
            end
            for k = 1:length(sigma_bac)
                state(j+1, 13+k) = state(j, 13+k) + normrnd(0,...
                    sqrt(At)*sigma_bac(k), [1, 1]);
            end
        end
    end

    % compute measurements
    for j = 1:Nt
         % extract bias
         bias = state(j, 8:16);

        [ddU_N] = compute_nBody(state(j, 1:3)' ,t(j), C_mat, S_mat, ...
                planetParams, posE(:, j), posM(:, j), posS(:, j));

        % rotate to body frame
        maxInd = 3 * j; minInd = maxInd -2;
        BN = BN_matrix(minInd:maxInd, :);
        ddU = BN * ddU_N * BN';

        T(1, j) = ddU(1,1) + bias(1) + noise_GG(1, j);
        T(2, j) = ddU(1,2) + bias(2) + noise_GG(2, j);
        T(3, j) = ddU(1,3) + bias(3) + noise_GG(3, j);
        T(4, j) = ddU(2,2) + bias(4) + noise_GG(4, j);
        T(5, j) = ddU(2,3) + bias(5) + noise_GG(5, j);
        T(6, j) = ddU(3,3) + bias(6) + noise_GG(6, j);

% %         T(1, j) = bias(1) - (0.9E-9)/measDim_QGG + noise_GG(1, j);

        r3 = state(j, 1:3)' - posS(:, j);
        eta= state(j, 7); 

         % compute SRP acceleration
        [aSRP, ~, ~] = SRP(r3, eta, mass, A,...
            planetParams);

        acc(:, j) = BN * aSRP + bias(7:9)' + noise_A(:, j);
    end
end

