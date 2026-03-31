clear;
clc;
close all;

%% COMPUTE RMS OF ERROR IN TIME DOMAIN FROM PSD VALUE

% sampling frequency [Hz]
fs = 1;   

% frequency range [Hz]
N = 1E-4;
f = 1E-5:N:fs/2;

% Instrument PSD [units^2 / Hz]
S_ASC = ((8.5E-6) * sqrt(f.^-1)).^2;                            % [rad^2 / Hz]
S_GRACE_FO = (2E-6)^2 * ones(1, length(f));                     % [rad^2 / Hz]
S_MICRO = (1E-10.*sqrt(0.4 + 0.001.*(f.^-1) + 2400.*(f.^4))).^2;  % [rad^2/s^2/Hz]
S_ASTRIX = (3E-8 * sqrt(1+4.6E-8 * (f.^-2))).^2;     % [rad^2/s^2/Hz]

% Select the sensor PSD
S = S_ASTRIX;

% plot 
figure()
loglog(f, sqrt(S)); xlabel('Hz'); ylabel('rad / \sqrt(Hz)'); grid on;
title('Instrument sqrt(PSD)');

% numerical integration
var   = trapz(f, S);   % numerical integration
sigma = sqrt(var);       % RMS noise
