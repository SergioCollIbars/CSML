clear;
clc;
close all;

addpath('functions/')
%%              GENERATE 1/F + WHITE NOISE
% Description: Generate noise profile combining the 1/f + white noise
% profile. Then, save theese values in a matlab file


% White noise sigma values;
sigma = [1.42E-12, 4.51E-11, 6.455E-11];               % [1/s^2]
measDim = 7.1040e-12;                                  % [1/s^2]

% Inputs
tmax = 1.1787e+06;       % [sec]
N    = 19645;            % Number of points time series

% Time and freq domains
x = linspace(1, tmax, N);
At = x(2) - x(1);
Fs = 1/At;
nfft = 2^nextpow2(N);
frec = 0 : Fs/(nfft-1) : Fs/2;

% noise vector
noise  = ones(3, N) * NaN;
noise_PSD = ones(3, length(frec)) * NaN;
for j = 1:length(sigma)
    sigma_wn = sigma(j);
    s = 6.32E-12 / sigma_wn;      % scaling value
    sigma_pn = sigma_wn / (50 * s);

    % generate noise
    MB = [5E-5, 0.1];
    [X, PSD] = noise_profile(sigma_wn, sigma_pn, N, Fs, MB);

    % save noise and PSD
    X = X./measDim;
    noise(j, :)  = X;
    noise_PSD(j, :) = PSD;
end

% plot combined noise
figure()
subplot(2, 1, 1)
plot(x./86400, noise)
title('Combined noise time series')
xlabel('Time [days]')

legend('\Gamma_{xx}', '\Gamma_{xy}', '\Gamma_{yz}')

subplot(2, 1, 2)
loglog(frec, sqrt(noise_PSD)/1E-9)
grid on;
title('Pink + White noise PSD')
xlabel('Frequency [Hz]')
ylabel('PSD^{1/2} [E /√Hz]')

% save time series
name = "noiseSeed_flicker_f60_T1.mat";
save(name, "noise");