clear;
clc;
close all;

f = 1/30;
tmax = 300 * 86400;
N = round(tmax * f);
TIME = linspace(0, tmax, N);
dt = TIME(2) - TIME(1);

S =(1E-3 / sqrt(2))^2;                                  % [E^2 / Hz]
fc= 1E-3;
q = S * 4 * sin(pi*fc/f)^2 * f^2;                       % [E^2 / s]

% bias
b_RW = ones(1, N);                                         % [E]

sigma_n = sqrt(S) * sqrt(f);                               % [E]
noise = normrnd(0, sigma_n, [1, N]);
q     = 1E-12;                                          % [E^2 / s]
for j = 1:N-1
    w2 = normrnd(0, sqrt(q*dt), [1, 1]);
    b_RW(j + 1) = b_RW(j) + w2;
end
b_RW = b_RW + noise;

figure()
plot(1:N, b_RW./1E-9 + noise./1E-9, 'LineWidth', 2)
ylabel('Eotvos')
title('Random Walk bias + noise time domain')

% compute PSD of the signal residual
res = b_RW - ones(1, N);
nfft     = [];                      % [] lets pwelch choose (or specify e.g. 2048)
window   = hamming(round(length(res)/8)); % or other window
noverlap = round(length(window)/2); % 50% overlap typical

% --- PSD estimation ---
[PSD, f] = pwelch(res, window, noverlap, nfft, f);

% --- Plot ---
figure;
loglog(f, PSD, 'color', [1, 0.41, 0.16]);
grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [E^2/Hz]'); 
title('PSD for gradiometer noise');