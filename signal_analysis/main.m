clear;
clc;
close all;
format long g;
%%                      MAIN CODE

% Paths
addpath('functions/')

% Inputs
N = 77760;            % Number of points time series
sigma_wn = 2.5E-10;   % STD GOCE (wn: White noise)
sigma_pn = 1.5E-10;   % STD GOCE (pn: Pink noise)

% Time and freq domains
x = linspace(1, 777600, N);
At = x(2) - x(1);
Fs = 1/At;
nfft = 2^nextpow2(N);
frec = 0 : Fs/(nfft-1) : Fs/2;

% generate noise
[X, PSD] = noise_profile(sigma_wn, sigma_pn, N, Fs);


% Plot
% combined noise
figure()
subplot(2, 1, 1)
plot(x./86400, X)
title('Combined noise time series')
xlabel('Time [days]')

subplot(2, 1, 2)
loglog(frec, sqrt(PSD)/1E-9)
grid on;
title('Pink + White noise PSD')
xlabel('Frequency [Hz]')
ylabel('PSD^{1/2} [E /√Hz]')