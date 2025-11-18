clear;
clc;
format long g;
close all;

addpath("NRHO_navigation/data/")
addpath("NRHO_navigation/functions/")
addpath("NRHO_navigation/functions/solver")
addpath("NRHO_navigation/functions/measurements")
addpath("NRHO_navigation/functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

%% METRIC COMPARING PERFORMANCE OF EKF vs UKF
% Description: Compare the advantage of UKF in the CR3BP problem using a
% statistic metric.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial configuration
system = "CR3BP";    % options: 2BP, CR3BP, F2BP, FCR3BP, EPHEM
frame  = "inertial"; % options: inertial, autonomous, RTN

% time parameters
T_orb = 1.3817;                        % [rad]
tmin = 0;                              % [rad]
tmax = 1 * T_orb;                      % [rad]
frec = 1/300;                           % [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME] = ...
    load_universe(system, [tmin, tmax], frec);

% load initial conditions. inertial frame (baricenter centered)
X0 = load_initCond(system, planetParams, TIME);

% compute integration time
n = round((TIME(end)-TIME(1))*(frec/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);

% load initial position and velocity uncertainty
sigmaP = 100E3;                                                   % [m]
sigmaV = 100;                                                     % [m/s]

sigmaP = sigmaP / planetParams(2);                              % [-]
sigmaV = sigmaV / (planetParams(2) * planetParams(3));          % [-]
P = diag([sigmaP, sigmaP, sigmaP, sigmaV, sigmaV, sigmaV].^2); % [m, m/s]

% Instrument measurement standardt deviation
sigmaM = 1E-9 * sqrt(frec);                                     % [1/s^2]
sigmaM = sigmaM/(planetParams(3)^2);                             % [-]
R = diag([sigmaM, sigmaM, sigmaM, sigmaM, sigmaM, sigmaM].^2);

%% integrate trajectory
disp('Computing trajectory ...')
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
STM0 = reshape(eye(6,6), [36, 1]);

[t, state_inertial] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0, 0), TIME, ...
    [X0; STM0], options);
TIME = t;

pos = state_inertial(:, 1:3); vel = state_inertial(:, 4:6); % [-]

%% compute non-liniarity indicators
disp('Computing non-liniarity indicators ...')
[Ns, Nm, trcP, trcV] = compute_nonLin_markers(pos, vel, P, R, TIME,...
    planetParams, poleParams, Cmat_true, Smat_true, system);

% clear kernel
cspice_kclear

%% compute lambda coefficient, upper and lower bound of gamma
disp('Compute upper & lower bound ...')
P = 0.8;                                    % expected minimum probability
a = ones(1, length(TIME)) * NaN;
b = ones(1, length(TIME)) * NaN;
for k = 1:length(TIME)
    f = 0.5 * (Ns(k) + Nm(k));
    lambda = min([f;1]);
    
    % upper bound 
    b(k) = lambda / sqrt(P);

    % lower bound
    a(k) = lambda/(P * sqrt(P)) * (2*(1 + sqrt(1-P))*(sqrt(P)-1) + P);
end
%% plot bounds for gamma
figure()
t = TIME./planetParams(3)/ 86400;   % [days]
a = a.*100; 
b = b.*100;
hold on
% Filled band
fill([t', fliplr(t')], [a, fliplr(b)], [0 0.447 0.741], ...
     'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName','Bounds');
plot(t', a, 'Color',[0 0.447 0.741], 'LineWidth',1.2, 'DisplayName','Lower');
plot(t', b, 'Color',[0 0.447 0.741], 'LineWidth',1.2, 'DisplayName','Upper');

xlabel('Time [days]')
ylabel('[%]')
title('\gamma upper & lower bounds')
grid on;

figure()
plot(t, 0.5.*(Ns+Nm),'LineWidth', 2)
xlabel('Time [days]')
grid on;
title('1/2 (N_m + N_s)')
figure()
plot3(pos(:, 1), pos(:, 2), pos(:, 3), 'LineWidth', 2)
axis equal;
grid on;
title('3D trajectory');

% periapsis & apoapsis bound value
pb = [max(a), max(b)];
ab = [min(a), min(b)];
disp('\Gamma value at periapsis: ' + string(pb))
disp('\Gamma value at apoapsis: ' + string(ab))

figure()
subplot(2, 1, 1)
plot(t, trcP * planetParams(2) / 1E3, 'LineWidth', 2)
ylabel('[km]')
grid on;
title('Position')
subplot(2, 1, 2)
plot(t, trcV * planetParams(2) * planetParams(3), 'LineWidth', 2)
xlabel('Time [days]');
title('Velocity')
ylabel('[m/s]')
grid on;

%% FUNCTIONS

function [Ns, Nm, trcP, trcV] = compute_nonLin_markers(pos, vel, P, R, TIME,...
    planetParams, poleParams, Cmat, Smat, system)
    
    sigmaCov = sqrt(diag(P)); sigmaV = sigmaCov(4); sigmaP = sigmaCov(1);
    Nt = length(TIME); Np = 2*length(P) + 1; dT = 60 * planetParams(3); %[-]

    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    STM0 = reshape(eye(6,6), [36, 1]);

    trcP = ones(1, Nt) * NaN; trcV = ones(1, Nt) * NaN;

    Ns  = ones(1, Nt) * NaN; Nm = ones(1, Nt) * NaN;
    fprintf('       Progress:    0%%');  % Initial message
    for k = 1:Nt
        % create 2Ns values
        X0 = [pos(k, :)'; vel(k, :)'];
        [~, state_0] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
                    poleParams, Cmat, Smat, system, 0, {0,0}, 0, 0), ...
                    [TIME(k), TIME(k)+dT], [X0; STM0], options);
        
        % dyanmics
        PHI_0 = reshape(state_0(end, 7:end), [6,6]);

        % measurement
        X0_prop = state_0(end, 1:6);
        [H_0] = compute_pos_partial_matrix(planetParams, poleParams, ...
                X0_prop, Cmat, Smat, TIME(k) + dT, eye(3,3));

        % compute Instant Measure Information Limit (IMIL)
% %         B = H_0' * (R\H_0);
% %         l = min(eig(B));
% %         poslimit = (1/l)^(0.5);
% %         posUnc = min([poslimit, sigmaP]);
% %         P0 = diag([posUnc, posUnc, posUnc, sigmaV, sigmaV, sigmaV].^2);
        P0 = diag([sigmaP,sigmaP,sigmaP, sigmaV, sigmaV, sigmaV].^2);
        stdP = sqrt(diag(P0));

        trcP(k) = sqrt(3) * stdP(1);
        trcV(k) = sqrt(3) * stdP(4);

        % function TBD --> return Xi vectors (6 x N)
        [Xi] = sigma_points(X0,P0,1.225,0);

        % integrate trajectory and compute indexes
        Ns_array = ones(1, Np) * NaN;
        Nm_array = ones(1, Np) * NaN;
        for j = 1:Np
            % dynamics non-liniearity index
            [~, state_i] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
                    poleParams, Cmat, Smat, system, 0, {0,0}, 0, 0), ...
                    [TIME(k), TIME(k)+dT], [Xi(:, j); STM0], options);
            PHI_i = reshape(state_i(end, 7:end), [6,6]);
            Ns_array(j) = norm(PHI_i - PHI_0, 'fro') / norm(PHI_0, 'fro');

            % measuremenr non-liniarity index
            Xi_prop = state_i(end, 1:6);
            [H_i] = compute_pos_partial_matrix(planetParams, poleParams, ...
                Xi_prop, Cmat, Smat, TIME(k) + dT, eye(3,3));

            Nm_array(j) = norm(H_i - H_0, 'fro') / norm(H_0, 'fro');
        end

        % select the maximum value at current time
        Ns(k) = max(Ns_array); Nm(k) = max(Nm_array);

        % Update every ~5% (optional)
        if mod(k, round(Nt/20)) == 0 || k == Nt-2
            fprintf('\b\b\b\b%3d%%', round(100 * k / Nt));
        end
    end
    fprintf('\n');
end

function [pos_earth, pos_moon] = computePos_circular(t, planetParams)
    mu = planetParams(1);
    M  = t;

    % Earth position
    pos_earth = -[mu*cos(M);mu*sin(M);0];           % [-]

    % Moon position
    pos_moon  = -[(mu-1)*cos(M);(mu-1)*sin(M); 0];  % [-]
end

function [H] = compute_pos_partials_CR3BP(planetParams, poleParams, x,...
    C_mat, S_mat, posE, posM, t, BN)
    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [ACI]
        rneg = x - Ar./2;   % [ACI]

        [ddU_pos] = compute_gradiometer_measurements(rpos, posE, posM, ...
                        planetParams, poleParams, C_mat, S_mat, t);

        [ddU_neg] = compute_gradiometer_measurements(rneg, posE, posM, ...
                        planetParams, poleParams, C_mat, S_mat, t);

        Ht = ((BN * ddU_pos * BN') - (BN * ddU_neg * BN'))./...
            (vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end

function [Hpos_mat] = compute_pos_partial_matrix(planetParams, poleParams, ...
    state, C_mat, S_mat, TIME, BN_mat)
    Hpos_mat = ones(6,3,length(TIME)) * NaN;

     % S/C position
    pos = state(:, 1:3)';   % [-]

    for k = 1:length(TIME)
         maxInd = 3* k; minInd = maxInd - 2;
         BN = BN_mat(minInd:maxInd, :); 
         % compute position Earth and Moon
        [posE, posM] = computePos_circular(TIME(k), planetParams);

        [H] = compute_pos_partials_CR3BP(planetParams, poleParams, pos(:, k),...
            C_mat, S_mat, posE, posM, TIME(k), BN);
        Hpos_mat(:, :, k) = H;
    end
end

function [X] = sigma_points(xhat,P,alpha,kappa)
    n = numel(xhat); L = chol(P,'lower');
    lambda = alpha^2*(n+kappa) - n;
    gamma  = sqrt(n + lambda);
    X = [xhat, xhat + gamma*L, xhat - gamma*L];   % n x (2n+1)
end

function [ddU] = compute_gradiometer_measurements(state, posE, posM, ...
    planetParams, poleParams, C_mat, S_mat, t)
% %     % rotation from Earth-Moon planet to J2000
% %     i_EM = deg2rad(0);   % [rad]
% %     EM_N = rotationMatrix(0, 0, i_EM, [1, 1, 1]);
% % 
% %     % extract rotation parameters
% %     RA_E = poleParams(1);            % RA Earth [rad]
% %     DEC_E = poleParams(2);           % DEC Earth [rad]
% %     W0_E = poleParams(3);            % prime meridian Earth [rad]
% %     WDot_E = poleParams(4);          % ang. velocity Earth [rad/s]
% % 
% %     RA_M = poleParams(5);            % RA Moon [rad]
% %     DEC_M = poleParams(6);           % DEC Moon [rad]
% %     W0_M = poleParams(7);            % prime meridian Moon [rad]
% %     WDot_M = poleParams(8);          % ang. velocity Moon [rad/s]
% % 
% %     Wt_E = WDot_E * t / planetParams(3) + W0_E;
% %     Wt_M = WDot_M * t / planetParams(3) + W0_M;
% %     ACAF1_N = rotationMatrix(pi/2 + RA_E, pi/2 - DEC_E, Wt_E, [3, 1, 3]);
% %     ACAF2_N = rotationMatrix(pi/2 + RA_M, pi/2 - DEC_M, Wt_M, [3, 1, 3]);
% %     
% %     ACAF1_EM = ACAF1_N * EM_N';
% %     ACAF2_EM = ACAF2_N * EM_N';

    ACAF1_EM = eye(3,3); ACAF2_EM = eye(3,3);

    % extract planet paramters (non-dimensional units)
    GM_E = planetParams(8) / (planetParams(2)^3 * planetParams(3)^2);     % [-]
    GM_M = planetParams(9) / (planetParams(2)^3 * planetParams(3)^2);     % [-]
    Re_E = planetParams(4) / planetParams(2);                             % [-]
    Re_M = planetParams(5) / planetParams(2);                             % [-]
    
    n_max      = planetParams(6);
    normalized = planetParams(7);

    relE = state(1:3) - posE;
    relM = state(1:3) - posM;
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, ~, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                ACAF1_EM*relE, Re_E, GM_E, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, ~, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                ACAF2_EM*relM, Re_M, GM_M, ...
                                                normalized);
    
    % Gravity gradient
    ddU =  (ACAF1_EM' * ddU1 * ACAF1_EM) + (ACAF2_EM' * ddU2 * ACAF2_EM);
end
