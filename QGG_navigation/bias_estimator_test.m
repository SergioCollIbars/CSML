clear;
clc;
close all;

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
tmax = 1*T + tmin;                                   % [-]
frec = 1/30; f_time = 1/30;

% define system
system = "CR3BP"; % options: 2BP, CR3BP, F2BP
[planetParams, ~, ~, ~, ~, ~] = ...
    load_universe(system, [tmin, tmax], frec);

% normalization values
measDim = planetParams(3)^2;
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
sigmaMeas = [1, 1/sqrt(2)] * 1E-12;                  % [1/s^2]
sigmaMeas = sigmaMeas./measDim;                      % [-]
R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
sigmaMeas(2), sigmaMeas(1)].^2);                     % [-]

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


Ax = 0; Nx = 0;
for k = 2:Nt
    PHI_1 = reshape(STM(k-1, :), [Ns, Ns]); % from 0 to k-1
    PHI_2 = reshape(STM(k, :), [Ns, Ns]);   % from 0 to k
    PHI = PHI_2 * inv(PHI_1);               % from k-1 to k;

    % compute measurements and visibility matrix
    [~, Ht, ~] = compute_measurements(TIME(k), state(k, 1:Ns), planetParams, ...
         poleParams, C_mat, S_mat, 0, 0, 0, [], DOM, posE(:, k), ...
         posM(:, k), posS(:, k), system);
    
    % gamma function
     At = t(k) - t(k-1);
     Gamma = At * [At/2*I;I];
     PHI_inv = inv(PHI);
     Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);

     PHIbd = exp(-dt/tau);         
     Sbd   = q * tau / 2 * (1 - exp(-2*dt/tau));       
     F = diag_stack(PHI, PHIbd, PHIbd, PHIbd, PHIbd, PHIbd, PHIbd);
     Qy = diag_stack(Gamma * qs * Gamma', Sbd, Sbd, Sbd, Sbd, Sbd, Sbd);
     M = inv(F)' * Aiprev_plus * inv(F);
     Ai_min = M - M*Qy*inv(eye(Nx, Nx) + M * Qy)*M;
     
     % measurement partials
     Hb = eye(6,6);
     h = [Ht, Hb];

    Ai_plus = Ai_min + h' * (R0(1:Nm, 1:Nm) \ h);

    obs(k) = rank(Ai_plus);

    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
end