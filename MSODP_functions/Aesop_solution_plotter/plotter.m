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
Input_GEO_1 = "/Users/sergiocollibars/Documents/Aesop_sol/500km/GEO.YY_215400";
Input_GEO_2 = "/Users/sergiocollibars/Documents/Aesop_sol/500km/GEO.XX_XZ_YY_ZZ_215400";

Input_true  = "/Users/sergiocollibars/Documents/Aesop_sol/GIF48.2007.GEO";

n_max = 120;
RE    = 6378136.3; % [m]

% Read GEO files
[C_1, S_1, sigma1_C, sigma1_S] = read_GEO_file(Input_GEO_1, n_max);
[C_2, S_2, sigma2_C, sigma2_S] = read_GEO_file(Input_GEO_2, n_max);
[C, S, sigma_C, sigma_S]       = read_true_GEO_file(Input_true, 360);

% Compute coeff. error
C_err1 = C(1:n_max+1, 1:n_max+1) - C_1; S_err1 = S(1:n_max+1, 1:n_max+1) - S_1;
C_err2 = C(1:n_max+1, 1:n_max+1) - C_2; S_err2 = S(1:n_max+1, 1:n_max+1) - S_2;

% compute RMS (GIF48 GEO)
rms_deg_true = rms_per_degree(C(1:n_max+1, 1:n_max+1), ...
    S(1:n_max+1, 1:n_max+1));
rms_cov_true = rms_per_degree(sigma_C(1:n_max+1, 1:n_max+1), ...
    sigma_S(1:n_max+1, 1:n_max+1));

% Hydrology level (mm)
H_lvl    = [4E-1,1.1E-1,5E-2,3E-2,2E-2,1.6E-2, 1E-2];
H_degree = [0,20,40,60,80,100,120];

% compute RMS (Standard)
rms_cov_err1  = rms_per_degree(sigma1_C, sigma1_S);
rms_deg_err1  = rms_per_degree(C_err1, S_err1);

% compute RMS (NSM)
rms_deg_err2  = rms_per_degree(C_err2, S_err2);
rms_cov_err2  = rms_per_degree(sigma2_C, sigma2_S);

% compute Geoid error (Estimated coefficients)
[degree_error_mm_1,degree_sigma_mm_1] = compute_geoid_error(n_max,...
    C_1, S_1, C(1:n_max+1, 1:n_max+1), S(1:n_max+1, 1:n_max+1), ...
    sigma1_C, sigma1_S, RE);

[degree_error_mm_2,degree_sigma_mm_2] = compute_geoid_error(n_max,...
    C_2, S_2, C(1:n_max+1, 1:n_max+1), S(1:n_max+1, 1:n_max+1), ...
    sigma2_C, sigma2_S, RE);

%% Plots

% plot Geoid error in a map
gridResolution = 1; % [deg]
C_t = C; S_t = S;
C_t(3, 1) = C(3, 1) + 0.4841668548961195E-03;
C_t(5, 1) = C(5, 1) - 0.7903040733333333E-06;
C_t(7, 1) = C(7, 1) + 0.1687251001365147E-08;

% True geoid
tt     = "Reference geoid value in m";
limVec = [-80, 80]; % [m]
plot_geoid(C_t(1:121, 1:121), S_t(1:121, 1:121), ...
    RE, gridResolution, tt, limVec);

% file 1 geoid error
tt     = "geoid error file 1 in cm";
limVec = [-5, 5]; % [cm]
plot_geoid(C_err1, S_err1, RE * 100, gridResolution, ...
    tt, limVec);

% file 2 geoid error
tt     = "geoid error file 2 in cm";
limVec = [-5, 5]; % [cm]
plot_geoid(C_err2, S_err2, RE * 100, gridResolution, ...
    tt, limVec);

% Plot Geoid error
figure();
semilogy(0:n_max, degree_error_mm_1, 'LineWidth', 1.5, 'Color','g')
hold on;
semilogy(0:n_max, degree_sigma_mm_1, '--', 'LineWidth', 1.5, 'Color','g')
semilogy(0:n_max, degree_error_mm_2, 'LineWidth', 1.5, 'Color','b')
semilogy(0:n_max, degree_sigma_mm_2, '--', 'LineWidth', 1.5, 'Color','b')
semilogy(H_degree, H_lvl, 'Color','k', 'LineWidth', 1.5, 'LineStyle','--');
grid on;
xlabel('Spherical harmonic degree')
ylabel('Geoid error [mm]')
legend('Error w.r.t. reference Standard', '', ...
    'Error w.r.t. reference NSM', '', 'Hydrology & Ice')

% Plot RMS error
figure();
xvals = (1:n_max+1) - 1; grayC = [0.5 0.5 0.5];
semilogy(xvals, rms_deg_true, 'LineWidth', 2, 'Color', 'k'); hold on;
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

% Plot pyramid graphs
figure();
semilogy(xvals, rms_deg_err1./rms_deg_err2, 'LineWidth', 1.5, 'Color', ...
    'r', 'LineStyle','none', 'Marker','square'); grid on;
title('Cov. ratio standard / NSM in RMS value');

plot_SH_coeffs(C_err1,S_err1, n_max, "SH error standard LS", [1E-15 5E-11], 'log');
plot_SH_coeffs(C_err2,S_err2, n_max, "SH error NSM LS", [1E-15 5E-11], 'log');

plot_SH_coeffs(C_err1./C_err2,S_err1./S_err2, n_max,...
    "SH error ratio Standard/NSM", [1E-4 1E4], 'log');

plot_SH_coeffs(sigma1_C./sigma2_C,sigma1_S./sigma2_S, n_max,...
    "\sigma ratio Standard/NSM", [], 'linear');

%% Stats
% % disp('C coefficients STATS: ');
% % compute_stats(abs(C_err1./C_err2), "C_{nm} coefficients");
% % disp('S coefficients STATS: ');
% % compute_stats(abs(S_err1./S_err2), "S_{nm coefficients}");

% % disp('C coefficients STATS: ');
% % compute_stats(abs(sigma1_C./sigma2_C), "C_{nm} coefficients");
% % disp('S coefficients STATS: ');
% % compute_stats(abs(sigma1_S./sigma2_S), "S_{nm coefficients}");