clear;
clc;
close all;

clear;
clc;
close all;
format long g;
addpath('../../NSM/functions/')
addpath('../../QGG_gravEstim/src/')
addpath('GOCE_L1b_MatlabReaders_1.1/data/')
addpath('GOCE_L2b_MatlabReaders/data/')
set(0,'defaultAxesFontSize',16);

% Earth planet data
path = "HARMCOEFS_EARTH_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
path = "SIGMACOEFS_EARTH_1.txt";
[sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
GM = 3.986004418E14;
n_max  = 10;
normalized = 1;
W = 0;                      % Rotation ang. vel   [rad/s]
W0 = 0;                     % Initial asteroid longitude
RA = -pi/2;                 % Right Ascension     [rad]
DEC = pi/2;                 % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

saveData = 0;

% trajectory data
date = "16_Nov_2012";
load(date + "_L2position.mat");
t = positions(:, 1);
pos_L2 = positions(:, 2:end);
Nt     = length(t);
t = t(1:Nt);
pos_L2 = pos_L2(1:Nt, :);
t_new = linspace(0, Nt, Nt);
pos_L2 = interpolate3DVector(pos_L2, t, t_new);

% Sampling setup
fs = t_new(2) -t_new(1);                  % Hz (1 sample per second)

% generate gradiometer data
sigma = 5E-12;                                      % [1/s^2]
noise = normrnd(0, sigma, [9, Nt]); noise0 = zeros(9, 1);
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, pos_L2, ...
                noise, Cnm, Snm, eye(3,3));
signal = Y(1, :);

% Add 1/f noise
noisy_signal = addOneOverFNoise(signal(1, :), fs);

% Filter between GOCE MBW: 0.005 to 0.1 Hz
filtered_signal = applyBandpassFilter(noisy_signal, fs, 0.005, 0.1);

% Plot time-domain signals
figure;
subplot(3,1,1); plot(t, signal); title('Original Signal'); ylabel('Amplitude');
subplot(3,1,2); plot(t, noisy_signal); title('Signal with 1/f Noise'); ylabel('Amplitude');
subplot(3,1,3); plot(t, filtered_signal); title('Filtered Signal (GOCE MBW)'); ylabel('Amplitude');
xlabel('Time (s)');

% Plot PSDs
figure;

% PSD of signal with 1/f noise
[pxx_noise, f_noise] = pwelch(noisy_signal, [], [], [], fs);
subplot(2,1,1);
loglog(f_noise, sqrt(pxx_noise)); grid on;
title('PSD of Signal with 1/f Noise');
xlabel('Frequency (Hz)'); ylabel('Amplitude Spectral Density');

% PSD of filtered signal
[pxx_filt, f_filt] = pwelch(filtered_signal, [], [], [], fs);
subplot(2,1,2);
loglog(f_filt, sqrt(pxx_filt)); grid on;
title('PSD of Filtered Signal');
xlabel('Frequency (Hz)'); ylabel('Amplitude Spectral Density');


%% FUNCTIONS
function noisy_signal = addOneOverFNoise(signal, fs)
% Adds GOCE-style noise to the input signal:
% - White noise in 0.005–0.1 Hz band
% - 1/f² behavior outside that band
% fs: sampling frequency in Hz

    N = length(signal);
    f = (0:N-1) * (fs/N);  % Frequency vector
    f(1) = 1e-6;           % Avoid division by zero

    % GOCE noise shaping parameters
    N0 = 1;       % white noise floor (arbitrary units, will normalize later)
    fc = 0.005;   % low cutoff frequency (Hz)
    fh = 0.1;     % high cutoff frequency (Hz)
    alpha = 2;    % exponent below fc
    beta = 2;     % exponent above fh

    % GOCE-like spectral shape
    noise_shape = sqrt(1 + (fc ./ f).^alpha + (f ./ fh).^beta);

    % Create white noise
    white_noise = randn(1, N);

    % FFT of white noise
    noise_fft = fft(white_noise);

    % Shape noise in frequency domain
    shaped_fft = noise_fft .* noise_shape;

    % Convert back to time domain
    goce_noise = real(ifft(shaped_fft));

    % Normalize power to match input signal
    goce_noise = goce_noise * std(signal) / std(goce_noise);

    % Add to signal
    noisy_signal = signal + goce_noise;
end


function filtered_signal = applyBandpassFilter(signal, fs, f_low, f_high)
% Applies a bandpass filter between f_low and f_high Hz
% fs: sampling frequency
% f_low, f_high: filter bounds in Hz

    % Butterworth bandpass filter
    [b, a] = butter(4, [f_low, f_high] / (fs / 2), 'bandpass');
    filtered_signal = filtfilt(b, a, signal);
end

function V_interp = interpolate3DVector(V, t_old, t_new)
% Interpolates a 3×N vector over time to a new time vector
% V: 3×N matrix (each row is a component of the vector)
% t_old: 1×N or Nx1 time vector for the original data
% t_new: 1×M or Mx1 new time vector

    V_interp = zeros(length(t_new), 3);  % Preallocate

    for i = 1:3
        V_interp(:, i) = interp1(t_old, V(:, i), t_new, 'pchip');  % or 'linear', 'spline'
    end
end
