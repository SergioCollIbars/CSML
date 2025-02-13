clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);
%%      NOISE GENERATION

% path
addpath("functions/")

% Inputs
tmax = 777600*1;
N = tmax/10;

% Time and freq domains
x = linspace(1, tmax, N);
At = x(2) - x(1);
Fs = 1/At;
nfft = 2^nextpow2(N);
frec = 0 : Fs/(nfft-1) : Fs/2;
ind_frec = find((5E-3 < frec) & (frec < 5.001E-3));

% noise std
sigma_wn = 6.32E-13;   % PSD GOCE
d = 6.32E-12 / sigma_wn;
sigma_pn = sigma_wn / (50 * d);   % PSD GOCE
%sigma_wn = 1.1E-8; % x PSD GRACE FO
%sigma_pn = 1E-8;   % x PSD GRACE FO
%sigma_wn = 3E-10;  % y PSD GRACE FO
%sigma_pn = 3E-12;  % y PSD GRACE FO

% pink noise function
a = pinknoise(1, N);
[PSD_a] = compute_PSD(fft(a, nfft), Fs, nfft, N);

% generate white noise.
STD = sigma_wn;
MEAN = 0;
wn = STD.*randn(N,1) + MEAN;
FFT_wn = fft(wn, nfft);
[PSD_wn] = compute_PSD(FFT_wn, Fs, nfft, N);
disp(PSD_wn(ind_frec))

% generate pink noise
[pn] = noiseGenerator_FIR(N, sigma_pn, 1);
PSD_pn = compute_PSD(fft(pn, nfft), Fs, nfft, N);
disp(PSD_pn(ind_frec))

% sum of both noise
[X, PSD] = noise_profile(sigma_wn, sigma_pn, N, Fs);
C = compute_PSD(fft(X, nfft), Fs, nfft, N);

% poyfit
p = polyfit(x, X, 1);
y = polyval(p,x);
res = X - y';


%% PLOT

% compatre noises
figure()
subplot(2, 1, 1)
loglog(frec, PSD_pn)
hold on;
loglog(frec, 1./(frec) / 1E26, 'LineWidth', 1.5, 'LineStyle','--' ,...
    'Color', 'k')
title('pn PSD')
grid on;
subplot(2, 1, 2)
loglog(frec, PSD_a)
hold on;
loglog(frec, 1./(frec) / 9E1, 'LineWidth', 1.5, 'LineStyle','--' ,...
    'Color', 'k')
title('a PSD')
grid on;

% White Noise time series
figure()
plot(x, wn)
title('White noise time series')

% Pink noise time series
figure()
plot(x, pn)
title('Pink noise time series')

% Combined noise time series
figure()
plot(x./86400, X)
title('Combined noise time series')
txt = "b = " + string(p(1)) + ", b_1 = " + string(p(2));
text(2E5,3.5E-11,txt,'FontSize',14)
xlabel('days')

% plot residuals
figure()
plot(x, res)
title('Combined noise residual')

% White noise PSD
figure()
loglog(frec, sqrt(PSD_wn)/1E-9)
hold on;
% % loglog(frec, ones(1, length(frec))*mean(sqrt(PSD_wn)/1E-9), ...
% %     'LineWidth', 1.5, 'Color', 'k', 'LineStyle','--')
title('White noise PSD')
xlabel('Frequency [Hz]')
ylabel('PSD [E /√Hz]')

figure()
loglog(frec, sqrt(PSD_pn)/1E-9)
xlabel('Frequency [Hz]')
ylabel('PSD [E/√Hz]')
grid on;
title('Pink noise PSD')
hold on;
loglog(frec, 1./(frec.^0.5) / 1E-6, 'LineWidth', 1.5, 'LineStyle','--' ,...
    'Color', 'k')
legend('PSD', '1/f trend')

figure()
loglog(frec, sqrt(PSD_pn)/1E-9)
xlabel('Frequency [Hz]')
ylabel('PSD [E /√Hz]')
grid on;
title('Pink noise PSD')

% combined noises
figure()
loglog(frec, sqrt(PSD)/1E-9)
grid on;
title('Pink + White noise PSD')
xlabel('Frequency [Hz]')
ylabel('PSD [E /√Hz]')
ax = gca;
ax.FontSize = 12;

% combined noises
figure()
c2 = load("c.mat").c;
c = sqrt(PSD)/1E-9;
smoothC = smoothdata(c,'movmean',500);
smoothC2 = smoothdata(c2,'movmean',500);
loglog(frec, c)
hold all;
loglog(frec, c2)
loglog(frec, smoothC, '--', 'Color','k', 'LineWidth', 1.5)
loglog(frec, smoothC2, '--', 'Color','k', 'LineWidth', 1.5)
grid on;
title('Required noise Power Spectral Density')
xlabel('Frequency [Hz]')
ylabel('PSD [E /√Hz]')
ax = gca;
ax.FontSize = 12;
legend('PSD \Gamma_{ij}', 'PSD \Gamma_{ii}', ...
    'Moving mean. Window 500 pt')

%% FUNCTION
function [Pxx] = compute_PSD(FFT, Fs, nfft, N)
    Pxx2 = abs(FFT).^2 /Fs /N; 
    Pxx2 = reshape(Pxx2, [nfft, 1]);
    Pxx = [Pxx2(1); 2*Pxx2(2:(nfft/2))];
end

function [X] = noiseGenerator_FIR(npts, sigma2, alpha)
    nn = 2*npts;
    ha = alpha/2;
    Qd = sqrt(sigma2);

    hfa = zeros(1, nn);
    wfa = zeros(1, nn);
    hfa(1) = 1;
    wfa(1) = Qd * randn();

    for i = 2:npts
        hfa(i) = hfa(i-1) * (ha+(i-2))/(i-1);
        wfa(i) = Qd * randn();
    end
    % the rest of the sequence is 0 by definition
    
    % perform FFT
    hfa = fft(hfa, npts);
    wfa = fft(wfa, npts);

    % multiply both
    wfa = wfa.*hfa;

    % do inverse FTT
    wfa = ifft(wfa, npts);
    X = wfa./npts;
end

function [X, PSD] = noise_profile(sigma_wn, sigma_pn, npts, Fs)
    % Description: compute the noise profile using white and pink noise
    % inputs: sigma_wn: STD white noise
    %         sigma_pn: STD pink noise
    %         npts: number of points
    %         Fs: sampling frequency
    % outputs: X: noise time series
    %          PSD: power spectral density
    
    % number of points FFT
    nfft = 2^nextpow2(npts);
    frec = 0 : Fs/(nfft-1) : Fs/2;

    % Gaussian noise time series
    X_wn = sigma_wn.*randn(npts,1);

    % Pink noise time series. FIR method
    [X_pn] = noiseGenerator_FIR(npts, sigma_pn, 1);
    
    % Band pass frequencies
    MB = [5E-3, 0.1];
    %MB = [5E-5, 0.1];
    [BP_FFT_wn] = bandPass_filter(fft(X_wn, npts), frec, MB);
    MB = [0, 5E-3];
    %MB = [0 , 5E-5];
    [BP_FFT_pn] = bandPass_filter(fft(X_pn, npts), frec, MB);
    BP_FFT = BP_FFT_pn + BP_FFT_wn;
    
    % compute noise from FFT to time series
    X = ifft(BP_FFT, npts);
    X = real(X);

    % total noise PSD
    PSD = compute_PSD(fft(X, nfft), Fs, nfft, npts);   
end


function [BP_FFT] = bandPass_filter(FFT, freq, MB)
    % Decription: Bandpass the fourier signal for the values specified in 
    %  the MB vector.
    % Input: FFT: values of th Fourier Transform
    %        freq: frequency vector
    %        MB: 2x1 vector with the band passs values
    % Ouput: BP_FFT: values of the band pass Fourier Transform

    ind1 = find((MB(1) < freq) & (freq < MB(2)));
    ind2 = length(FFT) + 1 - ind1;
    
    % Keep values inside the MB
    BP_FFT = zeros(length(FFT), 1);
    BP_FFT(ind1) = FFT(ind1);
    BP_FFT(ind2) = FFT(ind2);
end