clear;
clc;
close all;

addpath("NRHO_navigation/data/");
%% PLOT EART and MOON gravity field in RMS value

%% data text files
Earth_file_coeff  = "HARMCOEFS_EARTH_GGM05c.txt";
Moon_file_coeff   = "HARMCOEFS_MOON_GRGM660PRIM.txt";

Earth_file_uncrt  = "SIGMACOEFS_EARTH_GGM05c.txt";
Moon_file_uncrt   = "SIGMACOEFS_MOON_GRGM660PRIM.txt";

%% Read and plot Earth gravity field
file              = readmatrix(Earth_file_coeff);
Earth_n_max       = file(1); R_E = file(2);
[Nc, Ns, ~]       = count_num_coeff(Earth_n_max);

file              = readmatrix(Earth_file_coeff);
SH_coeff          = file(4:end);

file              = readmatrix(Earth_file_uncrt);
SH_uncrt          = file(4:end);

% compute RMS value
ZZ = zeros(Earth_n_max+1, Earth_n_max+1);
[Earth_RMS_SH_coeff] = computeRMS_coeffErr(Earth_n_max, Nc, Ns, SH_coeff, ...
    ZZ, ZZ);
[Earth_RMS_SH_uncrt] = computeRMS_coeffErr(Earth_n_max, Nc, Ns, SH_uncrt, ...
    ZZ, ZZ);

%% Read and plot Moon gravity field
file              = readmatrix(Moon_file_coeff);
Moon_n_max        = file(1); R_M = file(2);
[Nc, Ns, ~]       = count_num_coeff(Moon_n_max);

file              = readmatrix(Moon_file_coeff);
SH_coeff          = file(4:end);

file              = readmatrix(Moon_file_uncrt);
SH_uncrt          = file(4:end);

% compute RMS value
ZZ = zeros(Moon_n_max+1, Moon_n_max+1);
[Moon_RMS_SH_coeff] = computeRMS_coeffErr(Moon_n_max, Nc, Ns, SH_coeff, ...
    ZZ, ZZ);
[Moon_RMS_SH_uncrt] = computeRMS_coeffErr(Moon_n_max, Nc, Ns, SH_uncrt, ...
    ZZ, ZZ);

%% Plot gravity fields
figure()
semilogy(2:Earth_n_max, Earth_RMS_SH_coeff(2:end), 'LineWidth', 2, 'Color', 'b'); 
grid on; hold on;
semilogy(2:Moon_n_max, Moon_RMS_SH_coeff(2:end), 'LineWidth', 2, 'Color', 'g'); 
semilogy(2:Earth_n_max, Earth_RMS_SH_uncrt(2:end), 'LineWidth', 2, 'Color', 'b', ...
    'LineStyle', '--');
semilogy(2:Moon_n_max, Moon_RMS_SH_uncrt(2:end), 'LineWidth', 2, 'Color', 'g', ...
    'LineStyle', '--');
xlabel('Degree'); title('RMS value GGM05C and GRGM660PRIM');
legend('RMS GGM05C', 'RMS GRGM660PRIM',...
    '\sigma GGM05C', '\sigma GRGM660PRIM'); 
