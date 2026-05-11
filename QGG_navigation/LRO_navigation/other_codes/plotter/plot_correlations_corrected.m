clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);
addpath('/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/other_codes/error_budget/');
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');

%% PLOT CORRELATIONS
% Author: Sergio Coll-Ibars
% Description: Plot correlation coefficients for the first Monte Carlo
% realization using the current Low Lunar Simulator data structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data location. Use the same structure as plot_state.m
folder = "/Users/sergiocollibars/OneDrive - UCB-O365/Kepler_codes/results/LRO_March_20_2015_12H/";
file   = "orbit_LRO_SIM2_dataOut.mat";

% Load output file. Expected variables include:
%   time, mC_struct, I_ALIG, NB_MOON_mat
load(folder + file);

% Select only the first Monte Carlo simulation
fields = fieldnames(mC_struct);
if isempty(fields)
    error('mC_struct does not contain any Monte Carlo simulations.');
end

data = mC_struct.(fields{1});
Pf   = data.Pf;
Xf   = data.Xf;
Nt   = length(time);

% Time to UTC for plotting
utc  = cspice_et2utc(time', 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

% Compute rotation matrices using the same attitude/frame option as plot_state.m
BN_mat = compute_orientation_SC(time, Xf(1:6,:)', I_ALIG, NB_MOON_mat);

% Preallocate correlation arrays
C_pos_bias_x = nan(Nt, 12);
C_pos_bias_y = nan(Nt, 12);
C_pos_bias_z = nan(Nt, 12);
C_bias_SF    = nan(Nt, 6, 6);
sigma_pos    = nan(Nt, 3);

% Labels
biasSFLabels = {'b_{xx}', 'b_{xy}', 'b_{xz}', 'b_{yy}', 'b_{yz}', 'b_{zz}', ...
                'SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}'};
SFLabels     = {'SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}'};
posLabels    = {'radial', 'tangential', 'normal'};
biasLabels   = {'b_{xx}', 'b_{xy}', 'b_{xz}', 'b_{yy}', 'b_{yz}', 'b_{zz}'};

for k = 1:Nt
    p = reshape(Pf(k, :), [18, 18]);

    maxInd = 3*k;
    minInd = maxInd - 2;
    BN = BN_mat(minInd:maxInd, :);

    % Rotate position/velocity covariance into the spacecraft/body frame,
    % while leaving the remaining states in their original ordering.
    A   = blkdiag(BN, BN, eye(6), eye(6));
    P_B = A * p * A';

    sigma = sqrt(diag(P_B));
    C = P_B ./ (sigma * sigma');

    % Position correlations with gradiometer bias and scale-factor states
    C_pos_bias_x(k, :) = C(1, 7:end);
    C_pos_bias_y(k, :) = C(2, 7:end);
    C_pos_bias_z(k, :) = C(3, 7:end);

    % Bias-state correlations with scale-factor states
    C_bias_SF(k, :, :) = C(7:12, 13:18);

    sigma_pos(k, :) = sigma(1:3);
end

figure()
plot(tUTC, C_pos_bias_x, 'LineWidth', 2); grid on;
xlabel('Time [UTC]'); title('Position correlation: radial direction');
legend(biasSFLabels, 'Location', 'best');

figure()
plot(tUTC, C_pos_bias_y, 'LineWidth', 2); grid on;
xlabel('Time [UTC]'); title('Position correlation: tangential direction');
legend(biasSFLabels, 'Location', 'best');

figure()
plot(tUTC, C_pos_bias_z, 'LineWidth', 2); grid on;
xlabel('Time [UTC]'); title('Position correlation: normal direction');
legend(biasSFLabels, 'Location', 'best');

figure()
plot(tUTC, 3.*sigma_pos, 'LineWidth', 2); grid on;
xlabel('Time [UTC]'); title('3\sigma position formal uncertainty');
legend(posLabels, 'Location', 'best');

for k = 1:6
    figure();
    plot(tUTC, squeeze(C_bias_SF(:, k, :)), 'LineWidth', 2); grid on;
    xlabel('Time [UTC]');
    title(['Bias correlation to scale factors: ', biasLabels{k}]);
    legend(SFLabels, 'Location', 'best');
end

% Clear SPICE kernels
cspice_kclear
