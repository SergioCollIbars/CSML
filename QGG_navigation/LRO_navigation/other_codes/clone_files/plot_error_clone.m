clear;
clc;
close all;

%% PLOT SH ERROR IN CLONE FILES.

clone_folder     = "/Users/sergiocollibars/Desktop/CSML/matlab_codes/QGG_equations/clone_files/data_out/";
GRGM1200_SH_file = "/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data/HARMCOEFS_MOON_GRGM1200";
GRGM1200_SH_unct = "/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data/COEFSUNCRT_MOON_GRGM1200";

file             = readmatrix(GRGM1200_SH_file);
SH_coeff         = file(4:end);

file             = readmatrix(GRGM1200_SH_unct);
SH_uncrt         = file(4:end);

[Nc, Ns, Ncs]    = count_num_coeff(1200);

[Cnm, Snm]       = list2mat(1200, Nc, Ns, SH_coeff);

% get clone files
files  = dir(clone_folder);
files = files(~[files.isdir]);    % keep only files (remove folders)

n_max = 1200;
[RMS_SH_coeff] = computeRMS_coeffErr(n_max, Nc, Ns, SH_coeff, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_uncrt] = computeRMS_coeffErr(n_max, Nc, Ns, 1.*SH_uncrt, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

figure()
semilogy(1:n_max, RMS_SH_uncrt, 'LineWidth', 1.5, 'LineStyle', '--');hold on;
semilogy(1:n_max, RMS_SH_coeff, 'Color', 'k', 'LineWidth', 1.8); grid on;
for f = 1:length(files)
    clone  = readmatrix(clone_folder + string(files(f).name));
    clone  = clone(4:end);
    
   [RMS_SH_clone] = computeRMS_coeffErr(n_max, Nc, Ns, clone, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
   hold on;
   semilogy(1:n_max, RMS_SH_clone, '.', 'MarkerSize', 0.5);
end