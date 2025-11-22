%% Toy example: LS with GOCE-style G matrix
clear;clc;close all;
%% Parameters
N  = 4023;          % filter length and number of samples (odd)
dt = 1;             % sample interval [s]
fs = 1/dt;          % sampling frequency [Hz]
t  = (0:N-1)'*dt;

% Frequency grid (0 ... fs-df)
k  = (0:N-1)';
f  = k*fs/N;

%% Define Measurement Band (MB)
fL = 0.01;          % lower MB limit [Hz]
fH = 0.10;          % upper MB limit [Hz]
idxMB = (f>=fL & f<=fH);

%% 1/f + white noise PSD: 1/f below MB, flat inside MB
S_noise = ones(N,1);                % flat inside MB and above

idxLow = (f<fL & f>0);              % avoid f=0
S_noise(idxLow) = (fL./f(idxLow));  % 1/f below MB

% Enforce Hermitian symmetry for real time series
S_noise(N:-1:N/2+2) = S_noise(2:N/2);

%% Generate colored noise with this PSD
rng(1);
phi  = 2*pi*rand(N,1);                % random phase
Spec = sqrt(S_noise).*exp(1j*phi);
Spec(N/2+2:end) = conj(flipud(Spec(2:N/2)));  % Hermitian symmetry
noise = real(ifft(Spec))*sqrt(N);     % time-domain noise n(t)

%% Design zero-phase filter G as in the paper

% Desired "psd(fi)" of the filter: brick-wall in MB
psd_filt = zeros(N,1);
psd_filt(idxMB) = 1;

% Smooth transitions with a Hanning window of length 10
Lw = 10;
edges = find(diff(double(idxMB)) ~= 0);  % indices where band starts/ends

for e = 1:numel(edges)
    i0 = edges(e)-Lw+1;
    i1 = edges(e)+Lw;
    i0 = max(i0,1);
    i1 = min(i1,N);
    w  = hann(i1-i0+1);                % Hanning window

    if idxMB(edges(e)+1)               % transition 0 -> 1
        psd_filt(i0:i1) = max(psd_filt(i0:i1), w);
    else                               % transition 1 -> 0
        psd_filt(i0:i1) = min(psd_filt(i0:i1), w(end:-1:1));
    end
end

% Fourier matrix F (unitary DFT matrix)
n = (0:N-1)';
F = 1/sqrt(N) * exp(-1j*2*pi*(n*n.')/N);

% Filter matrix G
G = conj(F) * diag(sqrt(psd_filt)) * F;     % N x N

% (Optional) extract zero-phase FIR coefficients from central row
g = real(G((N+1)/2,:)).';

%% Check that G acts like a zero-phase convolution
noise_filt_mat = G*noise(:);                  % via matrix
noise_filt_conv = conv(noise, g, 'same');     % via convolution

% they should be almost identical
fprintf('max difference (matrix vs conv) = %.3e\n', ...
        max(abs(noise_filt_mat - noise_filt_conv)));

noise_filt = noise_filt_mat;   % use filtered noise

%% Simple LS problem: estimate amplitude of a sinusoid inside MB

f0 = 0.05;                      % signal frequency inside MB
H  = sin(2*pi*f0*t);            % N x 1 sensitivity
a_true = 1.0;

y = H*a_true + noise;          % raw measurements

% Apply G to both y and H (as in the paper)
y_tilde = real(G*y);
H_tilde = real(G*H);

% Ordinary LS on filtered system
a_hat = (H_tilde.'*H_tilde) \ (H_tilde.'*y_tilde);
fprintf('True amplitude = %.3f, LS estimate after filtering = %.3f\n', ...
        a_true, a_hat);

%% Plot noise spectra before and after filtering
Nfft = 4096;
[PSD_in, f_psd]  = pwelch(noise,      [], [], Nfft, fs);
[PSD_out, f_psd_2]     = pwelch(noise_filt, [], [], Nfft, fs);

figure;
loglog(f_psd, PSD_in,  'DisplayName','Raw noise'); hold on;
loglog(f_psd_2, PSD_out, 'DisplayName','Filtered noise');
xline(f0,   '--', 'DisplayName','Signal f_0');
xline([fL fH], ':',  'DisplayName','MB limits', 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('PSD');
legend('Location','best');
grid on;
title('Noise PSD before and after applying G');
