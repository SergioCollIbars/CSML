clear;
clc;
close all;

%%      FOURIER ANALYSIS OF THE 2ND ORDER TENSOR SIGNAL

% Input
addpath("data/")
set(0,'defaultAxesFontSize',16);

% data
phaseB_data = readtable('accData.txt');
GM = load("GM_periodogram.mat").A;
mission = "GOCE";       % GOCE / GRACE_FO

B_time = phaseB_data.TIME;
At = B_time(2) - B_time(1);
Fs = 1/At;
nfft = 2^nextpow2(length(B_time));
f = 0 : Fs/(nfft-1) : Fs/2;

T_xx = detrend(phaseB_data.ad_xx);
T_xy = detrend(phaseB_data.ad_xy);
T_xz = detrend(phaseB_data.ad_xz);
T_yy = detrend(phaseB_data.ad_yy);
T_yz = detrend(phaseB_data.ad_yz);
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

% signal noise
noise_ii = zeros(1, length(f));
noise_ij = zeros(1, length(f));
if(mission == "GRACE_FO")
    for j = 1:length(f)
        noise_ii(j) =  sqrt(1 + (f(j)/0.5)^4 + (0.1/f(j)));
        noise_ij(j) = 0.1 * sqrt(1 + (f(j)/0.5)^4 + (0.005/f(j)));
        if(f(j) < 5E-5)
             noise_ii(j)  = 1./f(j)/4E2;
             noise_ij(j)  = 1./f(j)/2E4;
        end
    end
elseif(mission == "GOCE")
    for j = 1:length(f)
        noise_ii(j) = 0.02;
        noise_ij(j) = 0.8;
        if(f(j) < 5E-3)
             noise_ii(j)  = sqrt(1./f(j)) / 7E2;
             noise_ij(j)  = sqrt(1./f(j)) / 1.8E1;
        end
    end
end
MB = linspace(1E-2, 1E4, length(f));

% plot noise 
figure()
loglog(f, noise_ii, 'LineWidth', 1.5, 'Color','k');
hold on;
loglog(f, noise_ij, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle','--');
xlabel('frequency [Hz]');
ylabel('noise [1/s^2/HZ]');
legend('x axis accelerometers', 'y and z axis accelerometers')
title('GRACE FO accelerometers noise profile');

% plot periodogram
figure();
subplot(2, 1, 1)
loglog(f, sqrt(P_xx)/1E-9, f, sqrt(P_yy)/1E-9, 'LineWidth', 1.5)
hold all;
plot(f, noise_ii, '--', 'color', 'k', 'LineWidth', 1.5)
if(mission == "GRACE_FO")
    plot(ones(1,length(f))*5E-5, MB, 'r--', 'LineWidth', 1.5)
elseif(mission == "GOCE")
    plot(ones(1,length(f))*5E-3, MB, 'r--', 'LineWidth', 1.5)
    plot(ones(1,length(f))*0.1, MB, 'r--', 'LineWidth', 1.5)
end
xlabel('frequency [Hz]')
ylabel('PSD^{1/2} [E / √HZ]')
title('Diagonal components')
legend('\Gamma_{xx}', '\Gamma_{yy}')

subplot(2, 1, 2)
loglog(f, sqrt(P_xy)/1E-9, f, sqrt(P_xz)/1E-9, f, sqrt(P_yz)/1E-9, ...
     'LineWidth', 1.5)
hold all;
plot(f, noise_ij, '--', 'color', 'k', 'LineWidth', 1.5)
if(mission == "GRACE_FO")
    plot(ones(1,length(f))*5E-5, MB, 'r--', 'LineWidth', 1.5)
elseif(mission == "GOCE")
    plot(ones(1,length(f))*5E-3, MB, 'r--', 'LineWidth', 1.5)
    plot(ones(1,length(f))*0.1, MB, 'r--', 'LineWidth', 1.5)
end
xlabel('frequency [Hz]')
ylabel('PSD^{1/2} [E / √HZ]')
title('Non diagonal components')
legend('\Gamma_{xy}', '\Gamma_{xz}', '\Gamma_{yz}')
sgtitle("Signal frequency analysis. Phase A")

% plot periodogram. All signals together
figure();
loglog(f, sqrt(P_xx)/1E-9, f, sqrt(P_yy)/1E-9, 'LineWidth', 1.5)
hold on;
loglog(f, sqrt(P_xy)/1E-9, f, sqrt(P_xz)/1E-9, f, sqrt(P_yz)/1E-9, ...
    'LineWidth', 1.5)
xlabel('frequency [Hz]')
ylabel('PSD [E / √HZ]')
title('Signal espectrum')
legend('\Gamma_{xx}', '\Gamma_{yy}', '\Gamma_{xy}', ...
    '\Gamma_{xz}', '\Gamma_{yz}')

% plot periodogram. All signals together without GM
figure();
loglog(f,abs(sqrt(P_xx) - sqrt(GM(:, 1)))/1E-9, ...
    f, abs(sqrt(P_yy) - sqrt(GM(:, 4)))/1E-9, 'LineWidth', 1.5)
hold on;
loglog(f, abs(sqrt(P_xy) - sqrt(GM(:, 2)))/1E-9, f, ...
    abs(sqrt(P_xz) - sqrt(GM(:, 3)))/1E-9, f, abs(sqrt(P_yz)-sqrt(GM(:, 5)))/1E-9, ...
    'LineWidth', 1.5)
xlabel('frequency [Hz]')
ylabel('PSD [E / √HZ]')
title('Signal espectrum without GM spectra')
legend('\Gamma_{xx}', '\Gamma_{yy}', '\Gamma_{xy}', ...
    '\Gamma_{xz}', '\Gamma_{yz}')

% plot time series of the signal
figure()
subplot(1, 2, 1)
plot(B_time./86400, T_xx, 'LineWidth', 1.5)
xlabel('Time [days]')
ylabel('1 / s^2')
title('\Gamma_{xx}')

subplot(1, 2, 2)
plot(B_time./86400, T_yy, 'LineWidth', 1.5)
xlabel('Time [days]')
ylabel('1 / s^2')
title('\Gamma_{yy}')

figure()
subplot(1, 3, 1)
plot(B_time./86400, T_xy, 'LineWidth', 1.5)
xlabel('Time [days]')
ylabel('1 / s^2')
title('\Gamma_{xy}')

subplot(1, 3, 2)
plot(B_time./86400, T_xz, 'LineWidth', 1.5)
xlabel('Time [days]')
ylabel('1 / s^2')
title('\Gamma_{xz}')

subplot(1, 3, 3)
plot(B_time./86400, T_yz, 'LineWidth', 1.5)
xlabel('Time [days]')
ylabel('1 / s^2')
title('\Gamma_{yz}')



%% FUNCTIONS
function [P] = compute_periodogram(data, Fs, N, nfft)
    Pxx2 = abs(data).^2 /Fs / N; 
    P = [Pxx2(1); 2*Pxx2(2:nfft/2)];
end