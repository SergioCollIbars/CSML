clear; clc; close all;
addpath(genpath('../functions/'))
addpath('../../QGG_gravEstim/src/')
addpath('../data/')
addpath('../../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/')
addpath('../../QGG_gravEstim/data_files/')

set(0,'defaultAxesFontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   GRACE ATTITUDE RESIDUAL NSM                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Use the GRACE trajectory and attitude to compute the
% gradiometer residuals based on the GM05c static Earth gravity model.
% Author: Sergio Coll-Ibars
% Date: 23/02/2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%                Mesurement mask 
%           xx xy xz yx yy yz zx zy zz
mask     =  [1, 1, 1, 0, 1, 1, 0, 0, 1]';
printMask_console(mask);

%% Load planet constants
[planetParams, poleParams, Kaula, r, Xtrue, Iner] = loadPlanet("Earth");
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 
[Cnm, Snm] = list2mat(n_max, Nc, Ns, Xtrue);

n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);

%% Load trajectory & planet orientation
disp('  Load S/C trajectory ...');
load('GRACE_traj.mat');

t = data(7, 1:end);          % time vector
Nt = length(t);
state_t = data(1:6, :)'; % trajectory vector

% compute planet orientation parameters
[ACAF_ACI] = compute_planetAtt(poleParams, t, 0, 0);

disp('  Computing attitude ... ')
[orientation, Mext] = SC_orientation(t, state_t, "ENU",...
    0, ACAF_ACI, Iner);

% attitude nominal value
theta    = orientation(1:3, :);  % Euler angle      [rad]
thetaDot = orientation(4:6, :);  % Euler angle rate [rad /s]

% compute Body to ACI rotation matrix
B_ACI_mat = ones(3*Nt, 3);
for j = 1:length(t)
    maxVal = 3 * j; minVal = maxVal -2;
    
    yaw  = theta(1, j); pitch = theta(2,j); roll = theta(3, j);
    B_ACI =rotationMatrix(yaw, pitch, roll, ...
        [3, 2, 1]);

    B_ACI_mat(minVal:maxVal, :) = B_ACI;
end

%% Compute state errors
disp('  Generating state errors ... ')
fs = 1 / (t(2) -t(1)); % Hz
plotCommand = 1;
[At1] = create_noise_from_PSD(fs, Nt/fs, plotCommand); % yaw error
[At2] = create_noise_from_PSD(fs, Nt/fs, plotCommand); % pitch error
[At3] = create_noise_from_PSD(fs, Nt/fs, plotCommand); % roll error

att_Err = [At1,At2,At3]';

angVel_nom = zeros(3, Nt); angVel_true = zeros(3, Nt);
angAcc_nom = zeros(3, Nt); angAcc_true = zeros(3, Nt);

%% Compute GG measurements (GRF frame)
disp('  Compure GG measurements ...');
[GG] = compute_GG_time_series_parallel(state_t, length(t), ACAF_ACI, Cnm,...
    Snm, n_max, Re, GM, normalized);

[Ytrue] = add_angularComponents(GG, theta, att_Err(1:3, :), angVel_nom,...
    angAcc_nom); 
[Ynom] = add_angularComponents(GG, theta, zeros(3, Nt), angVel_nom,...
    angAcc_nom);

dY  = Ytrue(logical(mask), :) - Ynom(logical(mask), :);

% plot gradiometer signal
if(length(Ytrue(:, 1)) > 6)
    plot_signal(Ytrue, t);
end

%% Apply NSM to eliminate 1st order orientation errors
disp('  Filter measurement residuals with NSM ...')
dY_NSM     = dY.*0;
Rot3 = reshape(B_ACI_mat.', 3, 3, []);  % size: 3 x 3 x Nt
parfor j = 1:Nt 
    maxVal = 3 * j; minVal = maxVal -2;
    B_ACI = Rot3(:, :, j)';
    
    [Hrot_grad] = compute_rotPartials_analy(GG(:, j), B_ACI);
    V_rot       = null(Hrot_grad(logical(mask), :)');
    P_rot       = V_rot * V_rot';

    dY_NSM(:, j)= P_rot * dY(:, j);
end

%% Plot GG time series residuals 
mask_name = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"];
figure(); tt = mask_name(logical(mask));
for k = 1:length(dY_NSM(:, 1))
    semilogy(t, abs(dY_NSM(k, :)./1E-9), "LineStyle", "none","Marker", ".", ...
        'DisplayName', tt(k)); hold on;
end
ylabel('Etovos'); grid on; title('Absolute val. NSM residuals'); legend();

%% Power Spectral Density of measurement residuals
Fs = 1 / (t(2) - t(1));                 % Hz (sampling rate)
pltt = "PSD GG residuals with NSM";
plotPSD_multiChannel(dY_NSM./1E-9, Fs, pltt, tt);

pltt = "PSD GG residuals without NSM";
plotPSD_multiChannel(dY./1E-9, Fs, pltt, tt);


%% FUNCTION
function [ddU_ACI] = compute_GG_time_series_parallel(state, Nt, RotPlanet, Cnm,...
    Snm, n_max, Re, GM, normalized)
    ddU_ACI = nan(9, Nt);
    x       = state(:, 1); y = state(:, 2); z = state(:, 3); 
    Rot3 = reshape(RotPlanet.', 3, 3, []);  % size: 3 x 3 x Nt
    parfor j = 1:Nt
        rt_ACI = [x(j);y(j);z(j)];
        
        % rotation matrix
        ACAF_ACI = Rot3(:, :, j)';
        rt_ACAF = ACAF_ACI * rt_ACI;

        [~, ~, T_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rt_ACAF, Re, GM, normalized);
        T_ACI = ACAF_ACI' * T_ACAF * ACAF_ACI;
        
        ddU_ACI(:, j) = [T_ACI(1,1); T_ACI(1,2) ; T_ACI(1,3); T_ACI(2,1);...
         T_ACI(2,2); T_ACI(2,3) ; T_ACI(3,1); T_ACI(3,2); T_ACI(3,3)];
    end
end