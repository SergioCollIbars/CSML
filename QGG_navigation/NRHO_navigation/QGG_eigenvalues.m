clear;
clc;
close all;
format long g;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")
%%                    QGG NAVIGATION EIGENVALUES
% Description: Compute the eigenvalues over time for the covariance matrix
% Author: Sergio Coll Ibars
% Date: 04/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
clear;
close all;
clc;

% specify simulation time, frequency and noise
tmin = 0;                                            % Phi = 0  [-]
Period = 1.4968;                                     % [-]
tmax = 4*Period + tmin;                              % [-]
frec = 1/30; f_time = 1/30;

% define system
system = "CR3BP"; % options: 2BP, CR3BP, F2BP
[planetParams, ~, ~, ~, ~, ~] = ...
    load_universe(system, [tmin, tmax], frec);

% normalization values
measDim = 1E-9;
timeDim = planetParams(3);                            % [1/s]
posDim  = planetParams(2);                            % [m]
velDim  = planetParams(3) * planetParams(2);          % [m/s]

% load parameters & initial conditions
[planetParams, poleParams, C_mat, S_mat, TIME, ~] = ...
    load_universe(system, [tmin, tmax], frec);
X0 = load_initCond(system, planetParams);

% compute measurement time
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);
N = f_time / frec;
DOM = ones(1, N) * NaN;
for h = 1:n/N
    val = N * (h - 1) + 1;
    DOM(h) = TIME(val);
end
dt = TIME(2) - TIME(1);                              % [-]

% Measurement weights
% % measScale     = planetParams(3)^2 / 1E-9 * 1E-3;                            % from [-] to [E]
% % positionScale = planetParams(2) * 1E-3;                                     % from [-] to [km]
% % velocityScale = planetParams(2)*planetParams(3) * 1E-3;                     % from [-] to [km/s]

positionScale = 1;                                                         
velocityScale = 1E1;                                                        
measScale     = 1E6;

% % positionScale = 1;                                                         
% % velocityScale = 1;                                                        
% % measScale = 1;


% scaling matrix
T = diag([positionScale, positionScale, positionScale, ...
    velocityScale, velocityScale, velocityScale, ...
    measScale, measScale, measScale, measScale, measScale ,measScale]);
T_inv = inv(T);

sigmaMeas = [1, 1] * 1E-12 * sqrt(frec);             % [1/s^2]
sigmaMeas = sigmaMeas./(planetParams(3)^2);          % [-]
R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
sigmaMeas(2), sigmaMeas(1)].^2);                     % [-]

%  state process noise
sigmaQ_s = 1E-11/ (planetParams(2)*planetParams(3)^2);% [-]
qs       = diag([sigmaQ_s, sigmaQ_s, sigmaQ_s].^2).*1;
I        = eye(3, 3);

% A priori uncertainty
Nx = 12; Ns = 6; Nm = 6;
P0 = diag(repelem(1E3, Nx));

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, C_mat, S_mat, system, 0, {0,0}, 0, 0), TIME, ...
    [X0; STM0], options);

STM  = state(:, Ns+1:Ns+Ns*Ns);
Nt  = length(TIME);

% primaries motion
[posE, posM, posS] = compute_posPrimaries(TIME, planetParams, system);

% time and Bias noise parameters
At         = t(2) - t(1);
At_sec     = At / planetParams(3);
q_E        = 1E-6;                                                            % [E / sqrt(sec)]
q          = (q_E * 1E-9 /(planetParams(3)^2 * sqrt(planetParams(3))))^2;     % [ - / sec]
Qb         = q * At * eye(6); 

% compute PCRB
PCRB       = ones(Nt, Nx*Nx) * NaN; P = PCRB;
sigmaP     = ones(Nx, Nt) * NaN;
PCRB(1, :) = zeros(1, Nx*Nx);
obs        = ones(1, Nt) * NaN;
eigenVals  = ones(Nx, Nt) * NaN;
condNum    = ones(1, Nt) * NaN;
for k = 2:Nt
    PHI_1 = reshape(STM(k-1, :), [Ns, Ns]); % from 0 to k-1
    PHI_2 = reshape(STM(k, :), [Ns, Ns]);   % from 0 to k
    PHI = PHI_2/(PHI_1);                    % from k-1 to k;

    % compute measurements and visibility matrix
    [~, Ht, ~] = compute_measurements(TIME(k), state(k, 1:Ns), planetParams, ...
         poleParams, C_mat, S_mat, 0, 0, 0, [], DOM, posE(:, k), ...
         posM(:, k), posS(:, k), system);
    Hb = eye(6);

    % scale measurement partials
    h = [Ht, Hb] * T;
    
    % Process noise
     Qrv = At.*[1/3*At^2*qs, 1/2*At*qs;...
        1/2*At*qs, qs];

     % previous information
     Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);


     PHIbd   = eye(6);
     F       = diag_stack(PHI, PHIbd);
     F       = T_inv * F * T; 
     Qy      = diag_stack(Qrv, Qb);
     Qy      = T_inv * Qy * T_inv';

     M       = inv(F)' * Aiprev_plus * inv(F);
     Ai_min  = M - M*Qy*inv(eye(Nx, Nx) + M * Qy)*M;


    Ai_plus = Ai_min + h' * (R0(1:Nm, 1:Nm) \ h);

    % compute eigenvalues
    if(rank(Ai_plus) == Nx)
         p = inv(Ai_plus); 
         Cov = (p + p.')/2;
         eigenVals(:, k) = eig(Cov);
         condNum(k) = cond(Cov);
    end

    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
end

% select last orbit period
time_T = t/Period;
idx = find(time_T > 3, 1);    % first index where t is greater than 3

% plot eigenvalues
figure()
for k = 1:Nx
    subplot(3, 4, k)
    semilogy(t(idx:end), eigenVals(k, idx:end), 'LineWidth', 2)
    grid on;
end

figure()
semilogy(t(idx:end), condNum(idx:end), 'LineWidth', 2)
grid on;


%% FUNCTIONS
function BigMatrix = diag_stack(varargin)
    % Determine total size
    total_size = sum(cellfun(@(M) size(M, 1), varargin));
    
    % Initialize the big matrix
    BigMatrix = zeros(total_size);
    
    % Fill the diagonal blocks
    idx = 1;
    for i = 1:length(varargin)
        M = varargin{i};
        s = size(M, 1);
        BigMatrix(idx:idx+s-1, idx:idx+s-1) = M;
        idx = idx + s;
    end
end

function [Hb] = compute_biasDrift_measPartials()
    Hb = zeros(6, 12);
    for j = 1:6
        maxInd = j * 2;
        minInd = maxInd - 1;
        Hb(j, minInd:maxInd) = [1, 0];
    end
end

function [PHIbd] = compute_STM_FOGM(t, tau)
    PHIbd = [1, tau*(1 - exp(-t/ tau));...
        0, exp(-t/tau)];
end

function [Sbd] = compute_PN_FOGM(t, tau, sigmaQ_b)
    S11 = tau^2 * ((1-exp(-2*t/tau)) - 4*(1-exp(-t/tau)) + 2*t/tau);
    S12 = tau * (1 - exp(-t/tau))^2;
    S21 = tau * (1 - exp(-t/tau))^2;
    S22 = (1 - exp(-2*t/tau));

    Sbd = sigmaQ_b * tau / 2 * [S11,S12;S21,S22];
end

