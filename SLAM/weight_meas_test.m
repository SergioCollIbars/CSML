clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);
%%                  TITLE WEIGHT MEASUREMENTS TEST
% Description: test to weight measurement such that they are make
% spectrally white.
% Author: Sergio Coll
% Date: 12/03/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% generate noise realization
tmax     = 3*3600;                % max time
fs       = 1/1;                   % sampling frequency [Hz]
Nt       = round(tmax * fs);       % number of points
f_min    = 5E-3;                   % cut frequnecy [Hz]
sigma_w  = 1E-3;                   % units / sqrt(Hz)

sigma_rw2       = 4 * sigma_w^2 * sin(pi * f_min / fs)^2;

stdVal          = sigma_w * sqrt(fs);
noise_white     = normrnd(0,stdVal, [1, Nt]);
eta             = sqrt(sigma_rw2) * randn(Nt,1);        % RW increments
noise_flicker   = cumsum(eta);                           % random walk

% total noise
noise = noise_flicker + noise_white';

% compute & plot PDS
window   = round(4096/1);
noverlap = window/2;
nfft     = window;

[PSD_rw, f] = pwelch(noise, window, noverlap, nfft, fs, "onesided");
figure; grid on; hold on;
set(gca, "XScale", "log", "YScale", "log")

loglog(f, PSD_rw, 'LineWidth', 1.4)

xlabel('Frequency [Hz]')
ylabel('PSD [units^2/Hz]')
title('PSD of Random Walk Noise');

% weight measurements
N = Nt;
R = zeros(N,N);
qb = sigma_rw2;
for k = 1:N
    for l = 1:N
        R(k,l) = (stdVal^2)*(k==l) + qb*min(k,l);
    end
end

L = chol(R, 'lower');
Wy = L \ noise;    % equivalent to inv(L)*y, but better conditioned

[PSD_Wy, f] = pwelch(Wy, window, noverlap, nfft, fs, "onesided");
figure; grid on; hold on;
set(gca, "XScale", "log", "YScale", "log")
hold on;
loglog(f, PSD_Wy, 'LineWidth', 1.4)

