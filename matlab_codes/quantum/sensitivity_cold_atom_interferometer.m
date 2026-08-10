clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);

%% SENSITIVITY COLD ATOM 
% Description: calculate and plot the sensitivity of several types of 
% cold-atom interferometers.
% Author: Sergio Coll-Ibars
% Date: 07/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate sensitivity

% Total propagation time [sec]
T  = 10;

% Sensitivity for QGG pathfinder [mE / sqrt(Hz)]
[ASD, ~] = sensitivity_QGG_pathfinder(T);
ASD_mE   = ASD / 1E-12;
disp('QGG pathfinder sensitivity at total propagation time of ' + string(T));
disp(ASD_mE);

% Sensitivity for QPI [mE / sqrt(Hz)]
[ASD, A] = sensitivity_QPI_gradiometer(T);
ASD_mE   = ASD / 1E-12;
disp('QPI gradiometer sensitivity at total propagation time of ' + string(T));
disp(ASD_mE);


%% Sensitivity scan

% total interrogation time scan values
T = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];

sigma_mE_QGG = nan(1, length(T));
ASD_mE_QGG   = nan(1, length(T));

sigma_mE_QPI = nan(1, length(T));
ASD_mE_QPI   = nan(1, length(T));
for k = 1:length(T)
    [ASD, sigma] = sensitivity_QGG_pathfinder(T(k));
    sigma_mE_QGG(k) = sigma / 1E-12;
    ASD_mE_QGG(k)   = ASD / 1E-12;

    [ASD, sigma] = sensitivity_QPI_gradiometer(T(k));
    sigma_mE_QPI(k) = sigma / 1E-12 / sqrt(100);
    ASD_mE_QPI(k)   = ASD / 1E-12;
end

figure()
plot(T, sigma_mE_QGG, 'LineWidth', 2, 'Color','b');grid on; hold on;
plot(T, sigma_mE_QPI, 'LineWidth', 2, 'Color','g');
xlabel('Sampling time [s]'); ylabel('mE', 'Interpreter','latex');
title('Instrument noise per interrogation time');
legend('QGG Pathfinder', 'QPI gradiometer');

figure()
plot(T, ASD_mE_QGG, 'LineWidth', 2, 'Color','b');grid on; hold on;
plot(T, ASD_mE_QPI, 'LineWidth', 2, 'Color','g');
xlabel('Sampling time [s]'); ylabel('mE / $\sqrt{Hz}$', 'Interpreter','latex');
title('Instrument ASD per interrogation time');
legend('QGG Pathfinder', 'QPI gradiometer');

%% AXULARY FUNCTIONS

function [ASD, sigma] = sensitivity_QGG_pathfinder(T)
    % Return the sensitivity for two  Mach-Zehnder atom interferometers.
    % https://link.springer.com/article/10.1140/epjqt/s40507-025-00338-1
    % Input: total interrogation time [s] T (split + recomb)
    % Assumption: T = T_cycle

    % Atom interferometer parameters
    C  = 0.5;                       % interferometer contrast
    L  = 0.3;                       % baseline separation [m]
    Nd = 1E5;                       % number of atoms
    ke = 8 * 2 * pi / (780 * 1E-9); % effective wavevector [1/m]

    % half of the total propagation time
    T_h   = T / 2;

    sigma = 2/C * 1/sqrt(Nd) * 1/(ke * L * T_h * T_h); % [s^-2]

    ASD   = sigma * sqrt(T);  % [s^(-2) * s^(0.5)]
end


function [ASD, sigma] = sensitivity_QPI_gradiometer(T)
    % Return the sensitivity for bloch-band interferometers.
    % QPI annual Review slides (005d_Pursuit_Precision_Performance.pdf)
    % Input: total interrogation time [s] T (split + recomb)
    % Assumption: T  = T_cycle
    
    % Rubidium wavelenght and mass
    lambda = 780.24 * 1E-9; % [m]
    mass   = 1.443 * 1E-25; % [km]
    
    % photom momentum
    h      = 6.62607015 * 1E-34; % [Js]
    h_bar  = h / (2*pi);
    hk     = h / lambda;

    % momentum 
    p = 8 * hk;

    % Fixed propagation times
    T1 = 100 * 1E-3; % [s]
    T2 = 200 * 1E-3; % [s]

    % Observable fringe [rad]
    delta_theta = pi / 2;
    
    % Number of effective atoms
    Neff = 10000;

    % sensitivity
    sigma = 1 / 8 * delta_theta * mass * h_bar * ...
            1 / (p^2 * T1 * T2 * (T2 + T)) * 1/sqrt(Neff); % [s^-2]

    ASD   = sigma * sqrt(T);  % [s^(-2) * s^(0.5)]
end