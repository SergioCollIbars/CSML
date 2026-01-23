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
%%                    QGG NAVIGATION OBSERVABILITY
% Description: compute the observability of the S/C state in the desired
% system using gradiometer measurements. Observability is function of the
% meas frequency and noise.
% Author: Sergio Coll Ibars
% Date: 04/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
clear;
close all;
clc;

% specify simulation time, frequency and noise
tmin = 0;                                            % Phi = 0  [-]
T = 1.4968;                                          % [-]
tmax = 4*T + tmin;                                   % [-]
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
sigmaMeas = [1, 1] * 1.5E-12 * sqrt(frec);             % [1/s^2]
sigmaMeas = sigmaMeas./measDim *1E-3 ;                     % [E]
R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
sigmaMeas(2), sigmaMeas(1)].^2);                     % [E^2]

%  state process noise
sigmaQ_s = (5E-11 / sqrt(3)) / (planetParams(2)*planetParams(3)^2); % [-]
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

% random-walk bias process noise
Nq         = 100;
sigmaQVec  = logspace(-7, -2, Nq).*1E-3;             % [E/ sqrt(sec)]
obsPercMat = ones(Nq, length(t));
RMS_val    = ones(3, Nq); max_val = ones(1, Nq); min_val = ones(1, Nq);
count      = 0;
At         = t(2) - t(1);
At_sec     = At / planetParams(3);
for j = 1:Nq
    q   = sigmaQVec(j)^2;                             % [E^2 / sec]

    % compute PCRB
    PCRB       = ones(Nt, Nx*Nx) * NaN; P = PCRB;
    sigmaP     = ones(Nx, Nt) * NaN;
    PCRB(1, :) = zeros(1, Nx*Nx);
    obs        = ones(1, Nt) * NaN;
    posUnc     = ones(3, Nt) * NaN;
    for k = 2:Nt
        PHI_1 = reshape(STM(k-1, :), [Ns, Ns]); % from 0 to k-1
        PHI_2 = reshape(STM(k, :), [Ns, Ns]);   % from 0 to k
        PHI = PHI_2/(PHI_1);                    % from k-1 to k;
    
        % compute measurements and visibility matrix
        [~, Ht, ~] = compute_measurements(TIME(k), state(k, 1:Ns), planetParams, ...
             poleParams, C_mat, S_mat, 0, 0, 0, [], DOM, posE(:, k), ...
             posM(:, k), posS(:, k), system);
        Ht = Ht * (planetParams(3)^2) / measDim *1E-3;   % [E/[-]]
        
         % Process noise
         Qrv = At.*[1/3*At^2*qs, 1/2*At*qs;...
            1/2*At*qs, qs];

         % previous information
         Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);
    
         PHIbd   = eye(6);
         Qb      = At_sec * q * eye(6);
         F       = diag_stack(PHI, PHIbd);
         Qy      = diag_stack(Qrv, Qb);
         M       = inv(F)' * Aiprev_plus * inv(F);
         Ai_min  = M - M*Qy*inv(eye(Nx, Nx) + M * Qy)*M;
         
         % measurement partials
         Hb = eye(6);
         h = [Ht, Hb];
    
        Ai_plus = Ai_min + h' * (R0(1:Nm, 1:Nm) \ h);
    
        obs(k) = rank(Ai_plus);
        if(obs(k)==12)
            p = inv(Ai_plus); sigmaP = sqrt(diag(p));
            posUnc(:, k) = sigmaP(1:3)*planetParams(2); % [m]
        end

        % store new value
        PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
    end

    % compute percentage of full observability
    obsPerc = obs;              % [-]

    obsPercMat(j, :) = obsPerc;
    count = count + 1;
    
    % select last orbit period
    time_T = t/T;
    idx = find(time_T > 3, 1);    % first index where t is greater than 3

    RMS_val(:, j) = rms(posUnc(:, idx:end), 2);
    normUnct      = vecnorm(posUnc(:, idx:end));
    max_val(:, j) = max(normUnct);
    min_val(:, j) = min(normUnct);

    disp('Progress = ' + string(count/(Nq) * 100) + ' %')
end

% compute error bars for the values
vals  = [1E-7 1E-6 1E-5 1E-4 1E-3];
errB  = zeros(2, length(vals));
yvals = zeros(1, length(vals)); 
for j = 1:length(vals)
    f = abs(sigmaQVec./1E-3 - vals(j));
    [~, idx] = min(f);
    
    errB(1, j) = max_val(idx);
    errB(2, j) = min_val(idx);

    yvals(j)   = vecnorm(RMS_val(:, idx));
end
figure()
loglog(sigmaQVec./1E-3, vecnorm(RMS_val), 'LineWidth', 2, 'Color', 'b');
hold all;
for k = 1:length(vals)
    line([vals(k) vals(k)], [errB(2, k) errB(1, k)], ...
        'LineWidth', 0.5, 'Color', 'b');
end
loglog(vals, yvals,    'ko', 'MarkerFaceColor','b', 'MarkerSize',6)   % nominal
loglog(vals, errB(1, :), 'ksq', 'MarkerFaceColor','b', 'MarkerSize',6)   % max
loglog(vals, errB(2, :), 'ksq', 'MarkerFaceColor','b', 'MarkerSize',6)   % min
grid on;
xlabel('bias rate [E $\sqrt{Hz}$]',  'Interpreter', 'latex');
ylabel('[m]');
xlim([1E-7 1E-3]); ylim([0, 10000]);
% % title('Formal error RMS value')
set(gca, 'YTick', [1 10, 100, 1000, 10000]);

% surface plot
figure()
time_h = t /planetParams(3) / 3600;
[X, Y] = meshgrid(time_h, sigmaQVec./1E-3); 
contourf(X, Y, obsPercMat, 'EdgeColor', 'none')
colormap("winter")                      % Specify colormap
c = colorbar;                                % Show color scale
c.Label.String = 'Number of states';   % <-- Your label here
xlabel('Time [hr]', 'Interpreter', 'latex')
ylabel('bias rate [E $\sqrt{Hz}$]',  'Interpreter', 'latex')
set(gca, 'YScale', 'log')              % Set Y axis to log scale

ax.GridColor = [0 0 0];   % darker grid
ax.GridAlpha = 1;         % more opaque
title('System observability over time');
view(2); xlim([1 90]);
xline(78.680000000000007, 'Color','r', 'LineStyle','--', 'LineWidth', 2);

% surface plot (v2)
figure(); clf
imagesc(time_h, sigmaQVec./1E-3, obsPercMat);  
axis xy                         
grid on
nStates = 12;
cmap = parula(nStates);   % or: lines(nStates)
colormap(cmap)

% Force discrete color mapping
clim([0.5 nStates + 0.5]);

cb = colorbar;
cb.Ticks = 1:nStates;
cb.TickLabels = string(1:nStates);
cb.Label.String = 'Number of observable states';

xlabel('Time [hr]', 'Interpreter', 'latex')
ylabel('bias rate [E $\sqrt{Hz}$]',  'Interpreter', 'latex')
set(gca, 'YScale', 'log')              % Set Y axis to log scale
xlim([1, 90]);
xline(78.680000000000007, 'Color','r', 'LineStyle','--', 'LineWidth', 2);
ylim([1E-6 1E-2]);
set(gca, 'YTick', [1E-6 1E-5, 1E-4, 1E-3, 1E-2])  % Specify 4 Y tick values

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

