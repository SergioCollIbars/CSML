clear; clc; close all;

%% Parameters
fs   = 1;           % [Hz] sampling frequency
N    = 5E4;        % number of samples (even)
t    = (0:N-1)'/fs; % time vector

% PSD model parameters
Swhite = 1E-1;   % white noise level [units^2/Hz]
fc     = 0.1;   % corner frequency between flat and 1/f region [Hz]
alpha  = 1;      % 1/f^alpha

%% 1) Build one-sided PSD model S_model(f)
df     = fs/N;
f_one  = (0:N/2)'*df;        % one-sided frequency grid [0 .. fs/2]

% Model: flat up to fc, then 1/f^alpha
f_min = 1e-5;                % avoid division by zero
f_eff = max(f_one, f_min);
S_one = Swhite^2 * (1 + (fc./f_eff).^alpha);  % [units^2/Hz]

% (Optional) plot the model PSD
figure;
loglog(f_one(2:end), S_one(2:end),'LineWidth',1.5);
xlabel('Frequency [Hz]');
ylabel('S_{model}(f) [units^2/Hz]');
title('Model PSD: white + 1/f');
grid on;

%% 2) Convert to two-sided PSD for length N (needed for FFT-based operations)
% S_one is one-sided. For a real signal:
%  - DC and Nyquist are single bins
%  - all interior positive freqs are mirrored -> two-sided is half there.
S_two_one              = S_one;
S_two_one(2:end-1)     = S_one(2:end-1)/2; % interior freqs halved
% Build full two-sided array (N points)
S_two_full = [S_two_one; flipud(S_two_one(2:end-1))];  % size N x 1

%% 3) Generate colored noise from the PSD model (frequency-domain shaping)

% Create frequency-domain coefficients consistent with S_two_full
% such that variance ~ integral S_one df. Then correct with a final scaling.
X = zeros(N,1);
% DC
X(1) = sqrt(S_two_full(1)*df)*randn;
% Nyquist (N/2+1)
X(N/2+1) = sqrt(S_two_full(N/2+1)*df)*randn;

% Positive freqs 2..N/2
for k = 2:N/2
    Sk    = S_two_full(k);
    sigma = sqrt(Sk*df/2);  % /2 because energy is shared with conjugate
    X(k)        = sigma*(randn + 1i*randn);
    X(N-k+1)    = conj(X(k));   % enforce Hermitian symmetry
end

% Time-domain noise realization
n = real(ifft(X));    % MATLAB ifft has 1/N internally

% Calibrate variance to match theoretical integral of S_one
var_theory = sum(S_one)*df;   % integral of one-sided PSD
var_emp    = var(n);
n          = n*sqrt(var_theory/var_emp);  % rescale

%% 4) Check that PSD of realization matches model (amplitude)
nfft_psd = N;
[PSD_est, f_est] = pwelch(n, hamming(nfft_psd/2), nfft_psd/4, nfft_psd, fs, 'onesided');

S_model_interp = interp1(f_one, S_one, f_est, 'linear', 'extrap');

figure;
loglog(f_est(2:end), PSD_est(2:end), 'LineWidth',1.5); hold on;
loglog(f_est(2:end), S_model_interp(2:end), '--', 'LineWidth',1.5);
xlabel('Frequency [Hz]');
ylabel('PSD [units^2/Hz]');
legend('Estimated PSD (pwelch)','Model PSD');
title('Comparison: model vs estimated PSD');
grid on;

%% 5) Build measurement model: y = H x_true + n
% Two parameters, e.g., offset and slow trend
x_true    = [1.0; -0.5];        % true parameters
H         = [ones(N,1), t/t(end)];  % N x 2

y         = H*x_true + n;       % measurements

%% 6) Unwhitened LS solution
x_LS = (H' * H) \ (H' * y);
C_LS = 3*sqrt(diag(inv(H'*H)));

%% 7) Build frequency-domain whitening "G matrix"
Gf = sqrt((Swhite) ./ S_two_full);   % frequency-domain whitening gain

%% 5) Apply zero-phase whitening to y and H
% Whitening in frequency domain:
Y  = fft(y);
Yw = Gf .* Y;                  % apply whitening gain
y_w = real(ifft(Yw));          % zero-phase whitening

Hw = zeros(size(H));
for j = 1:size(H,2)
    Hj  = fft(H(:,j));
    Hjw = Gf .* Hj;            % same filter applied to each column
    Hw(:,j) = real(ifft(Hjw));
end

% Whitened LS
x_WLS = (Hw' * Hw) \ (Hw' * y_w);
C_WLS = 3*sqrt(diag(inv(Hw'*Hw)));

%% 8) Compare solutions
disp('True parameters:');
disp(x_true.');

disp('Unwhitened LS estimate:');
disp(x_LS.');

disp('Whitened LS estimate (using G):');
disp(x_WLS.');

fprintf('Unwhitened error  : %g, %g\n', x_LS(1)-x_true(1), x_LS(2)-x_true(2));
fprintf('Whitened error    : %g, %g\n', x_WLS(1)-x_true(1), x_WLS(2)-x_true(2));
fprintf('Unwhitened 3 sigma  : %g, %g\n', C_LS(1), C_LS(2));
fprintf('Whitened 3 sigma    : %g, %g\n', C_WLS(1), C_WLS(2));

%% 9) (Optional) check that whitened noise is flat
NO = fft(n);
Nw = Gf.*NO;
n_w = real(ifft(Nw));

[PSD_w, f_w] = pwelch(n_w, hamming(nfft_psd/2), nfft_psd/4, nfft_psd, fs, 'onesided');

figure;
loglog(f_w(2:end), PSD_w(2:end),'LineWidth',1.5);
xlabel('Frequency [Hz]');
ylabel('PSD_{whitened} [units^2/Hz]');
title('PSD of whitened noise (should be ~flat)');
grid on;
