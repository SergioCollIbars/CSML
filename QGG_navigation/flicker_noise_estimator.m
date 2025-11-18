%% 1/f noise + RW-bias Kalman filter demo
clear; clc; close all;

%% Simulation parameters
fs   = 1;            % Hz (sampling frequency)
dt   = 1/fs;
N    = 2*5000;         % number of samples
t    = (0:N-1)'*dt;

alpha = 1;           % 1/f^alpha exponent (alpha = 1 => pink noise)
S0    = 1e-4;        % rough PSD scale factor for 1/f part (tune as needed)

sigma_v = 2e-3;      % measurement white noise std
q_rw    = 5.7e-6;      % RW process noise spectral density (tune!)

%% --- Step 1: Generate approximate 1/f noise in frequency domain ---
% Frequency grid
Nfft = N;
f    = (0:Nfft-1)'*(fs/Nfft);

% Avoid division by zero at DC: flatten PSD below f_min
f_min = fs/Nfft;
f_use = max(f, f_min);

% Target one-sided PSD ~ S0 / f^alpha over band
S_one = S0 ./ (f_use.^alpha);

% Build symmetric complex spectrum for real time series
% We'll populate bins 2..N/2 with random phase; enforce Hermitian symmetry.
X = zeros(Nfft,1);

% DC component (0 Hz): set to 0 mean for pure noise
X(1) = 0;

% Positive freqs (excluding Nyquist if Nfft even)
if mod(Nfft,2) == 0
    % Even length, Nyquist exists at bin Nfft/2+1
    k_pos = 2:(Nfft/2);      % strictly positive frequencies
    k_nyq = Nfft/2 + 1;
else
    % Odd length, no exact Nyquist bin
    k_pos = 2:((Nfft+1)/2);
    k_nyq = [];
end

% Random phases for positive freqs
phi = 2*pi*rand(numel(k_pos),1);
mag = sqrt(S_one(k_pos)/2);  % /2 because we will mirror spectrum

X(k_pos) = mag .* exp(1j*phi);

% Nyquist bin (if exists) must be real-valued
if ~isempty(k_nyq)
    X(k_nyq) = sqrt(S_one(k_nyq))*randn; 
end

% Negative freqs: Hermitian symmetry
if mod(Nfft,2) == 0
    % bins Nfft/2+2 ... Nfft mirror 2 ... Nfft/2
    X(Nfft/2+2:end) = conj(flipud(X(2:Nfft/2)));
else
    % bins (Nfft+3)/2 ... Nfft mirror 2 ... (Nfft+1)/2
    X((Nfft+3)/2:end) = conj(flipud(X(2:(Nfft+1)/2)));
end

% IFFT to time domain
b_1overf = real(ifft(X))*Nfft;   % scaling to keep power roughly consistent
b_1overf = b_1overf(:);
b_1overf = b_1overf - mean(b_1overf); % zero mean

% Optionally scale to a desired std
b_1overf = b_1overf * 1e-3;      % adjust overall amplitude

%% --- Step 2: Generate measurements (bias + white noise) ---
v   = sigma_v*randn(N,1);
z   = b_1overf + v;

%% --- Step 3: Kalman filter modeling bias as random walk ---
% Continuous-time RW q_rw [units^2/s] -> discrete Q = q_rw * dt
Q = q_rw * dt;
R = sigma_v^2;

x_hat = zeros(N,1);  % state estimate
P     = 1;           % initial variance (large-ish)

for k = 2:N
    % Time update
    x_pred = x_hat(k-1);     % x_{k|k-1} = x_{k-1|k-1}
    P_pred = P + Q;          % P_{k|k-1} = P_{k-1|k-1} + Q

    % Measurement update
    K      = P_pred / (P_pred + R);   % scalar Kalman gain
    x_hat(k) = x_pred + K*(z(k) - x_pred);
    P        = (1 - K)*P_pred;
end

%% --- Step 4: Evaluate performance ---
err = x_hat - b_1overf;
rmse = sqrt(mean(err.^2));

fprintf('RMSE of RW-KF estimate vs true 1/f bias: %.3e\n', rmse);

%% --- Step 5: Plots ---
figure;
subplot(3,1,1);
plot(t, b_1overf, 'LineWidth', 1); hold on;
plot(t, x_hat, '--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Bias');
legend('True 1/f bias','KF estimate');
title('1/f Bias and Random-Walk KF Estimate');

subplot(3,1,2);
plot(t, z, 'LineWidth', 0.8);
xlabel('Time [s]');
ylabel('Measurement');
title('Measurement (bias + white noise)');

subplot(3,1,3);
plot(t, err, 'k', 'LineWidth', 1); hold on;

% --- Compute covariance envelope (±3σ) ---
P_store = zeros(N,1);
P = 1;                 % initial variance
x_temp = 0;
P_store(1) = P;

for k = 2:N
    % Time update
    P_pred = P + Q;

    % Kalman gain
    K = P_pred / (P_pred + R);

    % Covariance update
    P = (1 - K) * P_pred;

    P_store(k) = P;       % store the updated covariance
end

sigma = sqrt(P_store);
plot(t,  3*sigma, 'r--', 'LineWidth', 1.2);
plot(t, -3*sigma, 'r--', 'LineWidth', 1.2);

xlabel('Time [s]');
ylabel('Error');
title('Estimation Error with ±3σ Covariance Envelope');
legend('Error','\pm3\sigma envelope');
grid on;
ylim([-4*sigma(100), 4*sigma(100)])


%% --- Optional: check PSD of true vs estimated bias ---
figure;
[PSD_true,f_w] = pwelch(b_1overf, hamming(512), 256, 1024, fs);
[PSD_est, f_w] = pwelch(x_hat,     hamming(512), 256, 1024, fs);

loglog(f_w, PSD_true, 'LineWidth', 1); hold on;
loglog(f_w, PSD_est, '--', 'LineWidth', 1);
xlabel('Frequency [Hz]');
ylabel('PSD');
legend('True 1/f bias','KF estimate');
title('PSD: True vs Estimated Bias');
grid on;


%% --- PSD of true 1/f, white noise, and combined measurement ---
figure;

% Compute PSDs
[PSD_1f,  f_w] = pwelch(b_1overf, hamming(512), 256, 1024, fs);
[PSD_w,   ~  ] = pwelch(v,         hamming(512), 256, 1024, fs);
[PSD_comb,~  ] = pwelch(z,         hamming(512), 256, 1024, fs);

% Plot
loglog(f_w, PSD_1f,  'LineWidth', 1.2); hold on;
loglog(f_w, PSD_w,   'LineWidth', 1.2);
loglog(f_w, PSD_comb,'LineWidth', 1.5);

xlabel('Frequency [Hz]');
ylabel('PSD');
title('PSD of 1/f Noise, White Noise, and Combined Signal');
legend('True 1/f noise','White noise','Combined (measurement)');
grid on;
