clear;
clc;
close all;

set(0, 'defaultAxesFontSize', 16);
rng(10);   % Reproducible noise realization

%% ================================================================
%  USER PARAMETERS
% ================================================================

% One-sided white-noise PSD [mE^2/Hz]
%
% Example:
% ASD = 1 mE/sqrt(Hz) = 1000 microE/sqrt(Hz)
ASD_w = 1;                 % [mE/sqrt(Hz)]
Sw    = ASD_w^2;             % [mE^2/Hz]

% White-noise / random-walk crossover frequency
fc = 1e-6;                   % [Hz]

% Measurement sampling/update frequency
fs = 1/300;                   % [Hz]

% Corresponding measurement cadence
dt = 1/fs;                   % [s]

% Simulation duration
Tsim = 2e5;                  % [s]

% Time vector
time = (0:dt:Tsim-dt)';
N    = length(time);

fprintf('White-noise ASD       = %.4g mE/sqrt(Hz)\n', ASD_w);
fprintf('White-noise PSD       = %.4g mE^2/Hz\n', Sw);
fprintf('Cutoff frequency      = %.4g Hz\n', fc);
fprintf('Sampling frequency    = %.4g Hz\n', fs);
fprintf('Sampling interval     = %.4g s\n', dt);
fprintf('Number of samples     = %d\n\n', N);

%% ================================================================
%  CHECK NYQUIST CONDITION
% ================================================================

if fs <= 2*fc
    warning(['The selected sampling frequency does not satisfy ' ...
             'the Nyquist condition fs > 2*fc.']);
end

%% ================================================================
%  RANDOM-WALK DRIVING PSD
%
%  Discrete random walk:
%
%       b_k = b_{k-1} + eta_k
%
%  Its PSD is:
%
%             S_RW
%  S_b(f) = -------------------------
%           4 sin^2(pi f / fs)
%
%  Choose S_RW so that S_b(fc) = Sw.
% ================================================================

SRW = 4 * Sw * sin(pi*fc/fs)^2;       % [mE^2/Hz]

fprintf('RW driving PSD        = %.4g mE^2/Hz\n', SRW);

%% ================================================================
%  CONVERT ONE-SIDED PSD TO DISCRETE-TIME VARIANCES
%
%  For a white sequence sampled at fs, with one-sided PSD S:
%
%       variance = integral_0^(fs/2) S df
%                = S*fs/2
% ================================================================

var_white = Sw  * fs/2;
var_RW    = SRW * fs/2;

sigma_white = sqrt(var_white);
sigma_RW    = sqrt(var_RW);

fprintf('White sample sigma    = %.4g mE\n', sigma_white);
fprintf('RW increment sigma    = %.4g mE/step\n', sigma_RW);
fprintf('sigma_RW/sigma_white  = %.4g\n\n', sigma_RW/sigma_white);

%% ================================================================
%  GENERATE NOISE REALIZATIONS
% ================================================================

% White measurement noise
white_noise = sigma_white * randn(N, 1);

% White driving sequence for the bias random walk
rw_increment = sigma_RW * randn(N, 1);

% Random-walk bias
bias_RW = cumsum(rw_increment);

% Total measurement noise
total_noise = white_noise + bias_RW;

%% ================================================================
%  TIME-DOMAIN PLOT
% ================================================================

Nplot = min(N, round(2e4/dt));

figure();

plot(time(1:Nplot)/3600, white_noise(1:Nplot), ...
    'LineWidth', 1.0);
hold on;

plot(time(1:Nplot)/3600, bias_RW(1:Nplot), ...
    'LineWidth', 1.5);

plot(time(1:Nplot)/3600, total_noise(1:Nplot), ...
    'LineWidth', 1.0);

grid on;
xlabel('Time [hr]');
ylabel('Noise [mE]');
title('White noise and random-walk bias realization');

legend('White measurement noise', ...
       'Random-walk bias', ...
       'Total noise', ...
       'Location', 'best');

%% ================================================================
%  ESTIMATE PSD USING WELCH METHOD
% ================================================================

segmentLength = min(4096, floor(N/4));
segmentLength = max(segmentLength, 256);

window   = hann(segmentLength);
noverlap = floor(segmentLength/2);
nfft     = max(4096, 2^nextpow2(segmentLength));

[Pwhite, f] = pwelch(white_noise, ...
    window, noverlap, nfft, fs, 'onesided');

[Prw, ~] = pwelch(bias_RW, ...
    window, noverlap, nfft, fs, 'onesided');

[Ptotal, ~] = pwelch(total_noise, ...
    window, noverlap, nfft, fs, 'onesided');

%% ================================================================
%  THEORETICAL PSD
% ================================================================

% Avoid the singular value at f = 0
f_theory = f;
f_theory(1) = NaN;

Sb_theory = SRW ./ ...
    (4*sin(pi*f_theory/fs).^2);

Stotal_theory = Sw + Sb_theory;

%% ================================================================
%  PSD PLOT
% ================================================================

figure();

loglog(f(2:end), Pwhite(2:end), ...
    'LineWidth', 1.1);
hold on;

loglog(f(2:end), Prw(2:end), ...
    'LineWidth', 1.1);

loglog(f(2:end), Ptotal(2:end), ...
    'LineWidth', 1.5);

loglog(f(2:end), Sw*ones(size(f(2:end))), ...
    '--', 'LineWidth', 1.8);

loglog(f(2:end), Sb_theory(2:end), ...
    '--', 'LineWidth', 1.8);

loglog(f(2:end), Stotal_theory(2:end), ...
    'k-', 'LineWidth', 2);

xline(fc, ':', 'f_c', ...
    'LineWidth', 1.8, ...
    'LabelVerticalAlignment', 'bottom');

grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [mE^2/Hz]');

title('White noise plus random-walk bias PSD');

legend('Simulated white noise', ...
       'Simulated RW bias', ...
       'Simulated total noise', ...
       'Theoretical white PSD', ...
       'Theoretical RW PSD', ...
       'Theoretical total PSD', ...
       'Location', 'best');

xlim([max(1/Tsim, f(2)), fs/2]);

%% ================================================================
%  ASD PLOT
%
%  This may be more convenient when comparing against an instrument
%  specification expressed in mE/sqrt(Hz).
% ================================================================

figure();

loglog(f(2:end), sqrt(Ptotal(2:end)), ...
    'LineWidth', 1.5);
hold on;

loglog(f(2:end), sqrt(Stotal_theory(2:end)), ...
    'k-', 'LineWidth', 2);

loglog(f(2:end), ASD_w*ones(size(f(2:end))), ...
    '--', 'LineWidth', 1.8);

xline(fc, ':', 'f_c', ...
    'LineWidth', 1.8, ...
    'LabelVerticalAlignment', 'bottom');

grid on;
xlabel('Frequency [Hz]');
ylabel('ASD [mE/\surdHz]');
title('Total-noise amplitude spectral density');

legend('Simulated total noise', ...
       'Theoretical total noise', ...
       'White-noise ASD', ...
       'Location', 'best');

xlim([max(1/Tsim, f(2)), fs/2]);

%% ================================================================
%  CONDITION: RW INCREMENT UNCERTAINTY < WHITE-NOISE UNCERTAINTY
%
%  sigma_RW / sigma_white
%
%       = sqrt(SRW / Sw)
%
%       = 2 |sin(pi fc/fs)|
%
%  Desired condition:
%
%       2 |sin(pi fc/fs)| < 1
% ================================================================

uncertainty_ratio = ...
    2 * abs(sin(pi*fc/fs));

condition_satisfied = uncertainty_ratio < 1;

fprintf('Condition 2|sin(pi fc/fs)| < 1:\n');

if condition_satisfied
    fprintf('SATISFIED for the selected fs.\n');
else
    fprintf('NOT SATISFIED for the selected fs.\n');
end

fprintf('Current uncertainty ratio = %.6f\n\n', ...
    uncertainty_ratio);

%% ================================================================
%  NUMERIC SEARCH FOR ACCEPTABLE fs VALUES
% ================================================================

fs_min = 2.001*fc;       % Just above Nyquist
fs_max = max(100*fc, 10*fs);

Nfs = 10000;

fs_test = logspace( ...
    log10(fs_min), log10(fs_max), Nfs);

fc_test = logspace( ...
    log10(fc), log10(2*fs), Nfs);

ratio_test = ...
    2 * abs(sin(pi*fc ./ fs_test));

ratio_test_2 = ...
    2 * abs(sin(pi*fc_test ./ fs));

valid_fs = (ratio_test < 1) & (fs_test > 2*fc);

%% Find the first acceptable fs
first_valid_index = find(valid_fs, 1, 'first');

if isempty(first_valid_index)
    fs_threshold_numeric = NaN;
    fprintf('No acceptable fs found in the selected search range.\n');
else
    fs_threshold_numeric = fs_test(first_valid_index);

    fprintf('Numerical minimum acceptable fs = %.8g Hz\n', ...
        fs_threshold_numeric);

    fprintf('Corresponding maximum cadence   = %.8g s\n', ...
        1/fs_threshold_numeric);
end

%% Analytical threshold in the Nyquist interval
fs_threshold_exact = 6*fc;

fprintf('Exact threshold fs > 6 fc       = %.8g Hz\n', ...
    fs_threshold_exact);

fprintf('Exact maximum cadence           = %.8g s\n\n', ...
    1/fs_threshold_exact);

%% ================================================================
%  PLOT GOOD fs REGION
% ================================================================

figure();

semilogx(fs_test, ratio_test, ...
    'LineWidth', 1.8);

hold on;

yline(1, '--', ...
    '\sigma_{RW}/\sigma_w = 1', ...
    'LineWidth', 1.8);

xline(2*fc, ':', ...
    'Nyquist: f_s = 2f_c', ...
    'LineWidth', 1.8);

xline(6*fc, '--', ...
    'Threshold: f_s = 6f_c', ...
    'LineWidth', 1.8);

xline(fs, '-.', ...
    'Selected f_s', ...
    'LineWidth', 1.8);

grid on;
xlabel('Sampling/averaging frequency, f_s [Hz]');
ylabel('\sigma_{RW}/\sigma_w');

title(['Condition for RW increment uncertainty ' ...
       'to remain below white noise']);

ylim([0, max(2, 1.05*max(ratio_test))]);

%% ================================================================
%  OPTIONAL: PLOT AS FUNCTION OF MEASUREMENT CADENCE
% ================================================================

Tsample_test = 1 ./ fs_test;

figure();

semilogx(Tsample_test, ratio_test, ...
    'LineWidth', 1.8);

hold on;

yline(1, '--', ...
    '\sigma_{RW}/\sigma_w = 1', ...
    'LineWidth', 1.8);

xline(1/(6*fc), '--', ...
    'T_s = 1/(6f_c)', ...
    'LineWidth', 1.8);

xline(dt, '-.', ...
    'Selected cadence', ...
    'LineWidth', 1.8);

set(gca, 'XDir', 'reverse');

grid on;
xlabel('Measurement cadence, T_s = 1/f_s [s]');
ylabel('\sigma_{RW}/\sigma_w');

title('Acceptable measurement cadence');

ylim([0, max(2, 1.05*max(ratio_test))]);

figure();
data = ratio_test_2;
idx  = fc_test > fs;
data(fc_test > fs) = NaN;
semilogx(fc_test, data, ...
    'LineWidth', 1.8);

grid on;
xlabel('Cut-off frequency value [Hz]');
ylabel('\sigma_{RW}/\sigma_w');
xlim([0, fs])
yline(1, '--', ...
    '\sigma_{RW}/\sigma_w = 1', ...
    'LineWidth', 1.8);