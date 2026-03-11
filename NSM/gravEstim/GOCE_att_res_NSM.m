clear; clc; close all;
addpath(genpath('../functions/'))
addpath('../../QGG_gravEstim/src/')
addpath('../data/')
addpath('../../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/')
addpath('../../QGG_gravEstim/data_files/')

set(0,'defaultAxesFontSize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   GOCE ATTITUDE RESIDUAL NSM                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Use the GOCE trajectory and attitude to compute the
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
load('Nov_L2position.mat');   % ECEF coordinates
load('Nov_L2velocity.mat');   % ECEF coordinates
load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF
load('Nov_L2ITRF2GRF.mat');   % rotation matrix ITRF 2 GRF

[~, t2]  = quaternion2CDM(outPut);
[~, t3]  = quaternion2CDM(outPut2);

% Check time-points to make sure that all files are at the same time.
t1 = positions(:, 1);
commonTimes = intersect(intersect(t1, t2), t3);

[~, idx1] = ismember(commonTimes, t1);
[~, idx2] = ismember(commonTimes, t2);
[~, idx3] = ismember(commonTimes, t3);

t0 = t1(idx1)';
t  = t0(1:end-1); 
Nt = length(t);

% % % WARNING: reducing data set
% % t = t(1:round(Nt/10));
% % Nt = length(t);

rn_ECEF  = positions(idx1, 2:end)';
vn_ECEF  = velocity(idx1, 2:end)';
ACI_ECEF = quaternion2CDM(outPut(idx2, :));
GRF_ACI  = quaternion2CDM(outPut2(idx3, :));

state_ACI = rotate2ECI([rn_ECEF; vn_ECEF], ACI_ECEF, t);
state_t(:, 1:3) = state_ACI(1:3, :)';
state_t(:, 4:6) = state_ACI(4:6, :)';

% compute planet orientation parameters
[ACAF_ACI] = compute_planetAtt(poleParams, t, 1, ACI_ECEF);

disp('  Computing attitude ... ')
[orientation, Mext] = SC_orientation(t, state_t, "GRF",...
    GRF_ACI, 0, Iner);

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
plotCommand = 0;
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
Rot2 = reshape(ACAF_ACI.', 3, 3, []);   % size: 3 x 3 x Nt
idx  = [4, 5, 6]; 

x       = state_t(:, 1); y = state_t(:, 2); z = state_t(:, 3); 
Ax_org  = 0; Ax_NSM = 0;
W       = eye(6) * (1E-12)^(-2);
ratio   = zeros(6, 3, Nt);
for j = 1:Nt 
    maxVal  = 3 * j; minVal = maxVal -2;
    B_ACI   = Rot3(:, :, j)';
    ACF_ACI = Rot2(:, :, j)';
    ACF_BODY = ACF_ACI * B_ACI';
    
    rt_ACI  = [x(j);y(j);z(j)];
    rt_ACAF = ACF_ACI * rt_ACI;

    [~, H_ACI]  = potentialGradient_Cnm(n_max, rt_ACAF, Re, ...
            GM, (ACF_ACI)', normalized);
    Hc_BODY     = rotate_coeffPartials(H_ACI, B_ACI);
    Hc          = Hc_BODY(logical(mask), 2:end);
    h_org       = Hc;

    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rt_ACI,...
        ACF_ACI, ACF_BODY);
    hpos   = Hpos(logical(mask), :);
  
    [Hrot_grad] = compute_rotPartials_analy(GG(:, j), B_ACI);
    [S, V, D]   = svd(Hrot_grad(logical(mask), :)');
    V_rot       = D(:, idx);
    P_rot       = V_rot * V_rot';

    hpos_nsm    = P_rot * hpos;
    r = hpos./hpos_nsm;
    % r1 = norm(r(:, 1)); r2 = norm(r(:, 2)); r3 = norm(r(:, 3));
    ratio(:, :, j) = r;

    dY_NSM(:, j)= P_rot * dY(:, j);
    h_nsm       = P_rot * Hc;

    Ax_org = Ax_org + h_org' * W * h_org;
    Ax_NSM = Ax_NSM + h_nsm' * W * h_nsm;
end

%% plot ratio between position error
figure()
tt = ["xx", "xy","xz", "yy", "yz", "zz"];
for k = 1:6
    subplot(2, 3, k)
    data = squeeze(ratio(k, :, :));
    semilogy(t, data, 'LineStyle','none', 'Marker','*'); grid on;
    title(tt(k));
    if(k==1), legend('Tangential', 'Normal', 'Radial'); end
end

%% plot information matrix ratio
r = diag(Ax_NSM)./diag(Ax_org);
figure()
plot(1:Ncs-1, r.*100, 'LineWidth', 2); grid on;
xlabel('Coefficient index'); title('Informaiton retention index');

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