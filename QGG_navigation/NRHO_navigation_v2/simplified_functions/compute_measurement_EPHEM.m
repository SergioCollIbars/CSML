function [T, B] = compute_measurement_EPHEM(t, state, planetParams,...
    BN_matrix, C_mat, S_mat, Nn, posE, posM, posS)
    % Compute measurements using the EPHEMERIDES model. The measurements
    % computed are: gradiometer measurements and acceloremeter measurements
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % S/C paramters
    mass = planetParams(11);    % [Kg]
    A    = planetParams(12);    % [-]
    eta  = 1.3;                 % [-]

    % Simulation points
    Nt  = length(t);
    
    % output matrices
    T   = ones(6, Nt) * NaN;
    acc = ones(3, Nt) * NaN;

    % Generate noise
    measDim_QGG = 1E-9;
    if(Nt > 1)
        dt = (t(2) - t(1)) / planetParams(3);                          % [se]
        frec = 1 / dt;                                                 % [Hz]
        At = (t(2) - t(1))/planetParams(3);                            % [sec]
    else
        frec = 1;
        At   = 1;
    end
    sigma_GG  = [1, 1, 1, 1, 1, 1] * 1E-3 * sqrt(frec) *1E-3;       % [killo-Eotvos (kE)]
    sigma_A   = 1E-8  * 1E3;                                        % [km/s^2]
    sigma_bGG = ones(1, 6) * 1E-6 * 1E-3 * sqrt(frec);              % [kE * sqrt(Hz)] 1E-8

    
    noise_GG = ones(6, Nt);
    for j=1:length(sigma_GG)
            noise_GG(j, :) = normrnd(0, sigma_GG(j), [1, Nt]) * Nn;
    end
    noise_A = normrnd(0, sigma_A, [3, Nt]) * Nn;

    % include bias random-walk
    if(Nn == 1)
        for j = 1:Nt-1
            for k = 1:length(sigma_bGG)
                state(j+1, 6+k) = state(j, 6+k) + ...
                    sqrt(At*(sigma_bGG(k)^2))*randn;
            end
        end
    end

    % compute measurements
    for j = 1:Nt
         % extract bias
         bias = state(j, 7:12);

        [ddU_N] = compute_nBody(state(j, 1:3)' ,t(j), C_mat, S_mat, ...
                planetParams, posE(:, j), posM(:, j), posS(:, j));

        % convert to 1/s^2
        ddU_N = ddU_N * (planetParams(3)^2);

        % rotate to body frame
        maxInd = 3 * j; minInd = maxInd -2;
        BN = BN_matrix(minInd:maxInd, :);
        ddU = BN * ddU_N * BN';

        T(1, j) = ddU(1,1)/measDim_QGG*1E-3 + bias(1) + noise_GG(1, j);
        T(2, j) = ddU(1,2)/measDim_QGG*1E-3 + bias(2) + noise_GG(2, j);
        T(3, j) = ddU(1,3)/measDim_QGG*1E-3 + bias(3) + noise_GG(3, j);
        T(4, j) = ddU(2,2)/measDim_QGG*1E-3 + bias(4) + noise_GG(4, j);
        T(5, j) = ddU(2,3)/measDim_QGG*1E-3 + bias(5) + noise_GG(5, j);
        T(6, j) = ddU(3,3)/measDim_QGG*1E-3 + bias(6) + noise_GG(6, j);

        r3 = state(j, 1:3)' - posS(:, j);

         % compute SRP acceleration
        [aSRP, ~, ~] = SRP(r3, eta, mass, A);

        acc(:, j) = BN * aSRP * 1E3  + noise_A(:, j);   % [km/s^2]
    end

    % return bias
    B = state(:, 7:12);
end

