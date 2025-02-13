function [T, Hmeas, Hcoeff] = compute_measurements(t, state, planetParams,...
    poleParams, C_mat, S_mat, Nn, Ncc, As, noiseSeed, DOM, posE, posM, posS,...
    system)
    %%                    COMPUTE MEASUREMENTS V2 FUNCTION
    % Description: compute position and velocity measurements without mask 
    % conditions given the truth state and time.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024
    % Data: Nn  = noise generation
    %       Ncc = Consider parameters 
    %       As  = Augmented state
    
    % Simulation points
    Nt = length(t);
    
    % output matrices
    T = ones(6, Nt) * NaN;
    
    % sensiticity matrix
    Hmeas = ones(6*Nt, 6) * NaN;
    
    % noise sigma
    measDim_QGG = planetParams(3)^2;
    measDim_Acc = planetParams(2)*planetParams(3)^2;
    sigma = [1, 1/sqrt(2), 1/sqrt(2), 1, 1/sqrt(2), 1] * 1E-12;  % [1/s^2]
    sigmaA = 1E-10;                                     % [m/s^2]
    sigma = sigma./measDim_QGG;                         % [-]
    sigmaA = sigmaA./measDim_Acc;                       % [-]

    % generate noise
    if(isempty(noiseSeed))
        flicker = 0;
        [noise]  = generate_noise(sigma, Nt, Nn, flicker, 0);
    else
        noise = noiseSeed * Nn;
    end
    
    % Augmented state
    if(As)
        T     = ones(9, Nt) * NaN;
        Hmeas = ones(9*Nt, 7) * NaN;

        noiseA = normrnd(0, sigmaA, [3, Nt]) * Nn;
        mass = planetParams(11);
        A    = planetParams(12);
    end

    % measurements computation
    for j = 1:Nt
        % compute gradiometer measurements. Inertial frame
        if(system == "EPHEM")
            [ddU] = compute_nBody(state(j, 1:3)' ,t(j), C_mat, S_mat, ...
                planetParams, posE(:, j), posM(:, j), posS(:, j));
            Hcoeff = 0;
        else
            [ACAF1_EM, ACAF2_EM] = compute_rotMat(poleParams, t(j));
            r1 = state(j, 1:3)' - posE(:, j);
            r2 = state(j, 1:3)' - posM(:, j);
            [ddU, Hcoeff] = compute_gradMeas(planetParams, ACAF1_EM, ACAF2_EM,...
                C_mat, S_mat, r1, r2, Ncc);
        end

        % cut off parameter
        % % Fn = sqrt(trace(ddU*ddU'));
        Fn = 0;
        m = isempty(find(t(j) == DOM, 1));

        if(Fn > 1E-9/measDim_QGG) || (m == 1)
            T(:, j) = NaN;
        else
            T(1, j) = ddU(1,1) + noise(1, j);
            T(2, j) = ddU(1,2) + noise(2, j);
            T(3, j) = ddU(1,3) + noise(3, j);
            T(4, j) = ddU(2,2) + noise(4, j);
            T(5, j) = ddU(2,3) + noise(5, j);
            T(6, j) = ddU(3,3) + noise(6, j);
        end

        % compute sensitivity matrix
        if(system == "EPHEM")
            [H] = compute_posPartials(planetParams, C_mat, S_mat, t(j),...
                state(j, 1:3)', posE(:, j), posM(:, j), posS(:, j));
        else
            [H] = compute_pos_par_CR3BP(planetParams, state(j, 1:3)', ...
                ACAF1_EM, ACAF2_EM, C_mat, S_mat, posE(:, j), posM(:, j));
        end

        if(As)
            r3 = state(j, 1:3)' - posS(:, j);
            eta= state(j, 7); 
            % compute SRP acceleration
            [aSRP, daSRP_dr, daSRP_deta] = SRP(r3, eta, mass, A,...
                planetParams);
            
            % measurement vector
            T(7:end, j) = aSRP + noiseA(:, j);
            if(m == 1), T(7:end, j) = NaN; end
            
            % store sensitivity matrix
            maxInd = 9*j;
            minInd = maxInd - 8;
    
            Hmeas(minInd:maxInd, :) = [H, zeros(6, 4);
                                        daSRP_dr, zeros(3, 3), daSRP_deta];
        else
            % store sensitivity matrix
            maxInd = 6*j;
            minInd = maxInd - 5;

            Hmeas(minInd:maxInd, :) = [H, zeros(6, 3)];
        end
    end
end

% compute gradiometer measurements
function [ddU, H] = compute_gradMeas(planetParams, ACAF1_EM, ACAF2_EM, ...
    C_mat, S_mat, r1, r2, Ncc)

    % extract planet paramters (non-dimensional units)
    GM1 = planetParams(8)./(planetParams(2)^3 * planetParams(3)^2);
    GM2 = planetParams(9)./(planetParams(2)^3 * planetParams(3)^2);
    Re1 = planetParams(4)./planetParams(2);
    Re2 = planetParams(5)./planetParams(2);

    n_max      = planetParams(6);
    normalized = planetParams(7);

    % compute gravity acceleration
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, ~, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                ACAF1_EM * r1, Re1, GM1, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, ~, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                ACAF2_EM * r2, Re2, GM2, ...
                                                normalized);

    % rotate back to inertial
    ddU1 = ACAF1_EM' * ddU1 * ACAF1_EM;
    ddU2 = ACAF2_EM' * ddU2 * ACAF2_EM;

    % total acceleration
    ddU = ddU1 + ddU2;
    
    % gravity tensor coefficient partials. Consider parameters
    if(Ncc)
        [~, h1] = potentialGradient_Cnm(n_max, ACAF1_EM*r1, Re1, GM1, ...
                ACAF1_EM', normalized);
        [~, h2] = potentialGradient_Cnm(n_max, ACAF2_EM*r2, Re2, GM2, ...
            ACAF2_EM', normalized);
        H = [h1(1:3, 2:end), h2(1:3, 2:end);...
             h1(5:6, 2:end), h2(5:6, 2:end); ...
             h1(9, 2:end), h2(9, 2:end)];
    else
        H = 0;
    end
end

function [ACAF1_EM, ACAF2_EM] = compute_rotMat(poleParams, et)
    % rotation from Earth-Moon planet to J2000
    i_EM = deg2rad(0);   % [rad]
    EM_N = rotationMatrix(0, 0, i_EM, [1, 1, 1]);

    RA_E = poleParams(1);            % RA Earth [rad]
    DEC_E = poleParams(2);           % DEC Earth [rad]
    W0_E = poleParams(3);            % prime meridian Earth [rad]
    WDot_E = poleParams(4);          % ang. velocity Earth [rad/s]

    RA_M = poleParams(5);            % RA Moon [rad]
    DEC_M = poleParams(6);           % DEC Moon [rad]
    W0_M = poleParams(7);            % prime meridian Moon [rad]
    WDot_M = poleParams(8);          % ang. velocity Moon [rad/s]

    Wt_E = WDot_E * et + W0_E;
    Wt_M = WDot_M * et + W0_M;
    ACAF1_N = rotationMatrix(pi/2 + RA_E, pi/2 - DEC_E, Wt_E, [3, 1, 3]);
    ACAF2_N = rotationMatrix(pi/2 + RA_M, pi/2 - DEC_M, Wt_M, [3, 1, 3]);

    ACAF1_EM = ACAF1_N * EM_N';
    ACAF2_EM = ACAF2_N * EM_N';
end

function [H] = compute_pos_par_CR3BP(planetParams, x, ACAF1_EM, ACAF2_EM, C_mat, S_mat, posE, posM)
    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [ACI]
        rneg = x - Ar./2;   % [ACI]

        r1 = rpos - posE; r2 = rpos - posM;
        [ddU_pos, ~] = compute_gradMeas(planetParams, ACAF1_EM, ACAF2_EM, ...
            C_mat, S_mat, r1, r2, 0);
        r1 = rneg - posE; r2 = rneg - posM;
        [ddU_neg, ~] = compute_gradMeas(planetParams, ACAF1_EM, ACAF2_EM, ...
            C_mat, S_mat, r1, r2, 0);
        Ht = (ddU_pos - ddU_neg)./(vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end

