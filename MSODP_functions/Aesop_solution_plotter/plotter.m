clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);
%% PLOT RESULTS FOR AESOP CODE
% Description: Given the file with actual error and nominal values, plot
% the results in RMS value and pyramid style.
% Author: Sergio Coll-Ibars
% Date: 04-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The values include the Scale Factor. (Amplitude scaled to 1E6)
Input_GEO_1 = "/Users/sergiocollibars/Documents/Aesop_sol/410km/GEO.XX_XZ_YY_ZZ_215000";
Input_GEO_2 = "/Users/sergiocollibars/Documents/Aesop_sol/410km/GEO.101101_0110_215400";

Input_true  = "/Users/sergiocollibars/Documents/Aesop_sol/GIF48.2007.GEO";

n_max = 120;

% Read GEO files
[C_1, S_1, sigma1_C, sigma1_S] = read_GEO_file(Input_GEO_1, n_max);
[C_2, S_2, sigma2_C, sigma2_S] = read_GEO_file(Input_GEO_2, n_max);
[C, S, sigma_C, sigma_S]       = read_true_GEO_file(Input_true, 360);

% % sigma2_C = (0.9).*sigma2_C; sigma2_S = (0.9).*sigma2_S;

% Compute coeff. error
C_err1 = C(1:n_max+1, 1:n_max+1) - C_1; S_err1 = S(1:n_max+1, 1:n_max+1) - S_1;
C_err2 = C(1:n_max+1, 1:n_max+1) - C_2; S_err2 = S(1:n_max+1, 1:n_max+1) - S_2;

% compute RMS (GIF48 GEO)
rms_deg_true = rms_per_degree(C(1:n_max+1, 1:n_max+1), ...
    S(1:n_max+1, 1:n_max+1));
rms_cov_true = rms_per_degree(sigma_C(1:n_max+1, 1:n_max+1), ...
    sigma_S(1:n_max+1, 1:n_max+1));

% compute RMS (Standard)
rms_cov_err1  = rms_per_degree(sigma1_C, sigma1_S);
rms_deg_err1  = rms_per_degree(C_err1, S_err1);

% compute RMS (NSM)
rms_deg_err2  = rms_per_degree(C_err2, S_err2);
rms_cov_err2  = rms_per_degree(sigma2_C, sigma2_S);

% plot
figure();
xvals = (1:n_max+1) - 1; grayC = [0.5 0.5 0.5];
semilogy(xvals, rms_deg_true, 'LineWidth', 2, 'Color', 'k'); hold all;
semilogy(xvals, rms_cov_true, 'LineWidth', 1.2, 'Color', grayC);
semilogy(xvals, rms_deg_err1, 'LineWidth', 1.5, 'Color', 'g', ...
    'Marker', '.', 'LineStyle','none'); grid on;
semilogy(xvals, rms_deg_err2, 'LineWidth', 1.5, 'Color', 'b', ...
    'Marker', '.', 'LineStyle','none');
semilogy(xvals, 3.*rms_cov_err1, 'LineWidth', 1.5, 'Color', 'g');
semilogy(xvals, 3.*rms_cov_err2, 'LineWidth', 1.5, 'Color', 'b');
title('RMS value AESOP solution');
legend('GIF48 GEO', '1\sigma ref. GEO', '', '', '3\sigma Standard LS', ...
    '3\sigma NSM LS');
% % legend('GIF48 GEO', '1\sigma ref. GEO', '', '', '2\sigma Att. error free', ...
% %     '2\sigma NSM LS');

figure();
semilogy(xvals, rms_cov_err1./rms_cov_err2, 'LineWidth', 1.5, 'Color', ...
    'r', 'LineStyle','none', 'Marker','square'); grid on;
title('Cov. ratio standard / NSM in RMS value');

plot_SH_coeffs(C_err1,S_err1, n_max, "SH error standard LS", [1E-15 5E-11], 'log');
plot_SH_coeffs(C_err2,S_err2, n_max, "SH error NSM LS", [1E-15 5E-11], 'log');

plot_SH_coeffs(C_err1./C_err2,S_err1./S_err2, n_max,...
    "SH error ratio Standard/NSM", [1E-4 1E4], 'log');

plot_SH_coeffs(sigma1_C./sigma2_C,sigma1_S./sigma2_S, n_max,...
    "\sigma ratio Standard/NSM", [], 'linear');

% stats
disp('C coefficients STATS: ');
compute_stats(abs(C_err1./C_err2), "C_{nm} coefficients");
disp('S coefficients STATS: ');
compute_stats(abs(S_err1./S_err2), "S_{nm coefficients}");

% % disp('C coefficients STATS: ');
% % compute_stats(abs(sigma1_C./sigma2_C), "C_{nm} coefficients");
% % disp('S coefficients STATS: ');
% % compute_stats(abs(sigma1_S./sigma2_S), "S_{nm coefficients}");