clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Compare Rummel's formulation vs the null space approach.
% Test in small radius orbits using Linear Leas Squares.
% Author: Sergio Coll
% Date: 11/10/25

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.36E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/10; dt = 1/f;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);
At = t(2) - t(1);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(5+Nc+Ns,5+Nc+Ns), [(5+Nc+Ns)^2, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
dX = ones(length(X)-1, 1) * 1E-1;

% Gravity estimation
sigma  = 1E-11 * sqrt(f);   % [10 mE / sqrt(Hz)]
R  =  diag([1,1,1,1,1,1])*sigma^2;
STD = sqrt(diag(R)');
mu = zeros(length(STD), 1);
Nm = length(R);
noise = normrnd(repmat(mu, 1, Nt), repmat(STD(:), 1, Nt));

% output variables
dY = zeros(6, Nt);
dY_1_NSM = zeros(3, Nt); dY_2_NSM = zeros(3, Nt); dY_3_NSM = zeros(3, Nt);
dY_1_PEP = zeros(6, Nt); dY_2_PEP = zeros(6, Nt); dY_3_PEP = zeros(6, Nt);
C_1_NSM  = zeros(3, 6, Nt); C_2_NSM  = zeros(1, 6, Nt); ...
C_3_NSM  = zeros(3, 6, Nt);
V_1_PEP  = zeros(6, 6, Nt); V_2_PEP  = zeros(6, 6, Nt); ...
V_3_PEP  = zeros(6, 6, Nt);
N = zeros(3, Nt);
for j = 1:Nt
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    rn_ACI = rn(:, j);
    
    % Inertially fixed
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

    w = [0;0;0];

    % computed meas. & partials
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
            zeros(9,1), Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(5, 2:end);...
         Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
    
     % compute Point mass approx error
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM,...
        rn_ACI, ACAF_ACI, ACAF_BODY);
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = ...
        compute_angularDyadPartials_v2(w, eye(3));
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];
    Hrot_grad = [Hrot_grad(1:3, :); Hrot_grad(5:6, :);Hrot_grad(9, :)];
    Hrot_omega_dyad = [Hrot_omega_dyad(1:3, :); ...
        Hrot_omega_dyad(5:6, :);Hrot_omega_dyad(9, :)];
    H_omegaDot_dyad = [H_omegaDot_dyad(1:3, :); ...
        H_omegaDot_dyad(5:6, :);H_omegaDot_dyad(9, :)];

    % nuisance parameter partials (3 cases)
    Hp_1 = Hpos;                                        % well conditioned
    Hp_2 = [Hpos, Hrot_grad];                           % ill-conditioned
    Hp_3 = Hrot_grad;                                   % well-conditioned

    % compute null space
    C_1 = null(Hp_1');
    [~,~,D] = svd(Hp_2');
    C_2 = D(:, 5);
    C_3 = null(Hp_3');
    
    % selected measurements
    Y = [Y_ACI(1:3);Y_ACI(5:6);Y_ACI(9)];
    dY(:, j) = Y; 

    % compute projected measurements
    dY_1_NSM(:, j) = (C_1' * Y);
    dY_2_NSM(:, j) = (C_2' * Y);
    dY_3_NSM(:, j) = (C_3' * Y);

    % compute projector
    PA1 = Hp_1 * pinv(Hp_1' * (R\Hp_1)) * (Hp_1'/(R));
    V_1 = eye(Nm,Nm) - PA1;
    PA2 = Hp_2 * pinv(Hp_2' * (R\Hp_2)) * (Hp_2'/(R));
    V_2 = eye(Nm,Nm) - PA2;
    PA3 = Hp_3 * pinv(Hp_3' * (R\Hp_3)) * (Hp_3'/(R));
    V_3 = (eye(Nm,Nm) - PA3);
    
    dY_1_PEP(:, j) = (V_1 * Y);
    dY_2_PEP(:, j) = (V_2 * Y);
    dY_3_PEP(:, j) = (V_3 * Y);

    % save null spaces
    C_1_NSM(:, :, j) = C_1';  C_2_NSM(:, :, j) = C_2';
    C_3_NSM(:, :, j) = C_3';

    V_1_PEP(:, :, j) = V_1;  V_2_PEP(:, :, j) = V_2;
    V_3_PEP(:, :, j) = V_3;

    % numerical noise
    N(:, j) = C_3'*noise(:, j);
end
% select test-case
testCase = 3;
if(testCase == 1)
    data_NSM = dY_1_NSM; data_PEP = dY_1_PEP;
    NSM = C_1_NSM; PEP = V_1_PEP;
elseif(testCase == 2)
    data_NSM = dY_2_NSM; data_PEP = dY_2_PEP;
    NSM = C_2_NSM; PEP = V_2_PEP;
elseif(testCase == 3)
    data_NSM = dY_3_NSM; data_PEP = dY_3_PEP;
    NSM = C_3_NSM; PEP = V_3_PEP;
end

% PSD parameters
win      = hamming(256);
noverlap = 128;             % 50%
nfft     = 512;

% compute PSD (Original)
[pxx, fs] = pwelch(dY(1, :), win, noverlap, nfft, f);
df = fs(2) - fs(1);

% compute noise
alpha = 2; fc = 3E-2;
[~, noise_f] = flicker_white(f, fs, sigma, fc, alpha);

pxx_mapped = zeros(length(pxx(:, 1)), length(dY(:, 1)));
Nxx_mapped = zeros(length(pxx(:, 1)), length(dY(:, 1)));
for ch = 1:length(dY(:, 1))
    [pxx_mapped(:, ch), ~] = pwelch(dY(ch, :),...
        win, noverlap, nfft, f);
    Nxx_mapped(:, ch) = noise_f;
end

% compute PSD (NSM)
[pxx, ~] = pwelch(data_NSM(1, :), win, noverlap, nfft, f);
pxx_mapped_NSM = zeros(length(pxx(:, 1)), length(data_NSM(:, 1)));
Nxx_mapped_NSM = zeros(length(pxx(:, 1)), length(data_NSM(:, 1)));
for ch = 1:length(data_NSM(:, 1))
    [pxx_mapped_NSM(:, ch), ~] = pwelch(data_NSM(ch, :),...
        win, noverlap, nfft, f);
    [A, ~] = pwelch(N(ch, :),...
        win, noverlap, nfft, f);

    for j = 1:6
        [val, ~] = pwelch(squeeze(NSM(ch, j, :)),...
            win, noverlap, nfft, f);
        [~, Z1] = psd_product_one_sided(val, noise_f, fs);
        Nxx_mapped_NSM(:, ch) =  Nxx_mapped_NSM(:, ch) + ...
            Z1./2;
    end
end

% compute PSD (PEP)
[pxx, ~] = pwelch(data_PEP(1, :), win, noverlap, nfft, f);
pxx_mapped_PEP = zeros(length(pxx(:, 1)), length(data_PEP(:, 1)));
Nxx_mapped_PEP = zeros(length(pxx(:, 1)), length(data_PEP(:, 1)));
for ch = 1:length(data_PEP(:, 1))
    [pxx_mapped_PEP(:, ch), ~] = pwelch(data_PEP(ch, :),...
        win, noverlap, nfft, f);
    for j = 1:6
        [val, ~] = pwelch(squeeze(PEP(ch, j, :)),...
            win, noverlap, nfft, f);
        [~, Z1] = psd_product_one_sided(val, noise_f, fs);
        Nxx_mapped_PEP(:, ch) =  Nxx_mapped_PEP(:, ch) + ...
            Z1./2;
    end
end

% plot noise
figure
loglog(fs, noise_f)
grid on;
title('PSD noise')

% Plot original signal vs projected signal ratio
figure;
semilogy(fs, sum(pxx_mapped_NSM, 2)./sum(pxx_mapped, 2), ...
    fs, sum(pxx_mapped_PEP, 2)./sum(pxx_mapped, 2), 'Linewidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');
title('Signal PSD ratio (Welch Method)');
legend('NSM', 'PEP')

% Plot original signal vs projected signal
figure()
semilogy(fs, sum(pxx_mapped, 2), fs, sum(pxx_mapped_NSM, 2), ...
    fs, sum(pxx_mapped_PEP, 2), 'Linewidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');
title('Signal PSD (Welch Method)');
legend('Original','NSM', 'PEP')

% Plot noise
figure()
loglog(fs, sum(Nxx_mapped_NSM, 2), ...
       fs, sum(Nxx_mapped_PEP, 2), 'Linewidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');
title('Noise PSD (Welch Method)');
legend('NSM', 'PEP')

% Plot SNR
figure()
semilogy(fs, sum(pxx_mapped, 2)./sum(Nxx_mapped, 2),...
    fs, sum(pxx_mapped_NSM, 2)./sum(Nxx_mapped_NSM, 2), ...
    fs, sum(pxx_mapped_PEP, 2)./sum(Nxx_mapped_PEP, 2), 'Linewidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');
title('SNR PSD (Welch Method)');
legend('Original','NSM', 'PEP')


%% functions
function [f, Z1] = psd_product_one_sided(A1, B1, f)
% A1, B1: one-sided PSDs on f = 0..fs/2 (column vectors, same size)
% f:      frequency vector (0..fs/2)
% fs:     sampling rate

df  = f(2)-f(1);
N1  = numel(f);                  % one-sided length (includes DC and Nyquist)
% Build even two-sided spectra (exclude duplicate DC & Nyquist)
A2  = [A1; A1(end-1:-1:2)];
B2  = [B1; B1(end-1:-1:2)];
N2  = numel(A2);                 % should be 2*N1 - 2

% Circular convolution over one DTFT period (keep PSD units)
Z2  = ifft( fft(A2).*fft(B2), 'symmetric') * (df);   % [N2 x 1]

% Fold to one-sided: DC and Nyquist not doubled, interior bins add mirrors
Z1        = zeros(N1,1);
Z1(1)     = Z2(1);               % DC
Z1(N1)    = Z2(N1);              % Nyquist
idx_pos   = 2:N1-1;              % interior positive freqs
idx_neg   = N2 - idx_pos + 2;    % their negative-frequency mirrors
Z1(idx_pos) = Z2(idx_pos) + Z2(idx_neg);
end

function [f, S] = flicker_white(fs, f, sigma_w, fc, alpha)
% FLICKER_WHITE  Simple 1/f^alpha + white PSD model (one-sided)
%
%   [f,S] = flicker_white(fs, Nf, sigma_w, fc, alpha)
%
% Inputs:
%   fs       - sampling frequency [Hz]
%   Nf       - number of frequency points (one-sided)
%   sigma_w  - standard deviation of white noise
%   fc       - corner frequency [Hz]
%   alpha    - slope exponent (1 for 1/f, 2 for 1/f^2, etc.)
%
% Output:
%   f        - frequency vector [Hz] (0..fs/2)
%   S        - one-sided PSD [units^2/Hz]
%
% Example:
%   [f,S] = flicker_white(0.1, 4097, 1e-9, 1e-3, 1);
%   loglog(f,S), grid on

    S0 = 2*sigma_w^2/fs;             % one-sided white level
    S = S0 * ones(size(f));          % start with flat spectrum

    % 1/f^alpha below fc
    idx = f <= fc & f > 0;
    S(idx) = S0 * (fc ./ f(idx)).^alpha;
    S(1) = S(find(f>0,1));           % avoid Inf at DC
end
