clear;
clc;
close all;

% Input
addpath("data/")
set(0,'defaultAxesFontSize',16);

% data
data = readmatrix('measData.txt');
TIME = data(:, 1);

xx = data(:, 2);
xy = data(:, 3);
xz = data(:, 4);
yy = data(:, 5);
yz = data(:, 6);

At = TIME(2) - TIME(1);
Fs = 1/At;
nfft = 2^nextpow2(length(TIME));
f = 0 : Fs/(nfft-1) : Fs/2;

T_xx = xx;
T_xy = xy;
T_xz = xz;
T_yy = yy;
T_yz = yz;
N = length(T_xx);

% fourier espectra
FFT_xx = fft(T_xx, nfft);
FFT_xy = fft(T_xy, nfft);
FFT_xz = fft(T_xz, nfft);
FFT_yy = fft(T_yy, nfft);
FFT_yz = fft(T_yz, nfft);

% compute periodogram
[P_xx] = compute_periodogram(FFT_xx, Fs, N, nfft);
[P_xy] = compute_periodogram(FFT_xy, Fs, N, nfft);
[P_xz] = compute_periodogram(FFT_xz, Fs, N, nfft);
[P_yy] = compute_periodogram(FFT_yy, Fs, N, nfft);
[P_yz] = compute_periodogram(FFT_yz, Fs, N, nfft);


%%              PLOTS

% plot periodogram
figure();
subplot(2, 1, 1)
loglog(f, P_xx, f, P_yy, 'LineWidth', 1.5);
legend('\Gamma_{xx}', '\Gamma_{yy}')
ylabel('PSD [-]')

subplot(2, 1, 2)
loglog(f,P_xy, f, P_xz, f, P_yz, 'LineWidth', 1.5)
xlabel('frequency [-]')
ylabel('PSD [-]')
sgtitle('Gradiometer, periodogram of the signal')
legend('\Gamma_{xy}', '\Gamma_{xz}', ...
    '\Gamma_{yz}')
%%              FUNCTION
function [P] = compute_periodogram(data, Fs, N, nfft)
    Pxx2 = abs(data).^2 /Fs / N; 
    P = [Pxx2(1); 2*Pxx2(2:nfft/2)];
end