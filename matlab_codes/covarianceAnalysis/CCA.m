clear;
clc;
close all;

%%              CONSIDER COVARIANCE ANALYSIS

% paths
addpath('functions/')
addpath('data/')

%  load planet parameters
[Cmat, Smat] = readCoeff('HARMCOEFS_BENNU_OSIRIS_0.txt');
Re = 246;           % [m]
GM = 5.2;           % [m^3 / s^2]

% load position and rotation matrices
filename = 'orbitData.txt';
posData = importdata(filename).data;
r_B = [posData(:, 5), posData(:, 6), posData(:, 7)]';
v_B = [posData(:, 11), posData(:, 12), posData(:, 13)]';

% read position in ACI coords
r_ACI = [posData(:, 2), posData(:, 3), posData(:, 4)]';

% Rotation matrices
filename = 'ECI2Body.txt';
BN = importdata(filename).data;

filename = 'N2ACAF.txt';
ACAF_N = importdata(filename).data;

% input time simulation parameters
tmin= 0;                % [s]
tmax = 777600;          % [s]
N = tmax / 10;          % [-]
time = linspace(tmin, tmax, N);

% input STD for state and consider parameters
sigmaX_ii = 6.32E-13;
sigmaX_ij = 6.32E-13;
sigmaC_b = 0;           % 0
sigmaC_bdot = 2E-15;    % 2E-15

% compute variances
sigmaX2_ii = sigmaX_ii*sigmaX_ii;
sigmaX2_ij = sigmaX_ij*sigmaX_ij;
sigmaC2_b = sigmaC_b*sigmaC_b;
sigmaC2_bdot = sigmaC_bdot*sigmaC_bdot;


Rx = diag([sigmaX2_ii, sigmaX2_ij, sigmaX2_ij, sigmaX2_ij, sigmaX2_ii, ...
        sigmaX2_ij, sigmaX2_ij, sigmaX2_ij, sigmaX2_ii]);
Rc =diag([sigmaC2_b, sigmaC2_bdot]); 


% START COVARIANCE ANALYSIS
Ax = 0;
Mx = 0;
for k = 1:length(time)
    maxInd= 3*k;
    minInd = maxInd - 2;
    R_ACAF_N = ACAF_N(minInd:maxInd, :);

    % current position in the ACAF frame
    r_ACAF = R_ACAF_N * r_ACI(:, k);

    % compute state sensitivity matrix.
    [H] = potentialGradient_Cnm(6, r_ACAF, Re, GM, ...
        R_ACAF_N', 0);
    
    % compute consider sensitivity matrix.
    hc = [1, time(k)];
    Hc = [hc;hc;hc;hc;hc;hc;hc;hc;hc];

    % sum covariances
    Ax = Ax + H' *  inv(Rx) * H;
    Mx = H' * inv(Rx) * Hc;
end
Px = inv(Ax);
Mx = -Px * Mx;
Pc = Px + Mx * Rc * Mx';

sigmaX = sqrt(diag(Px));
sigmaC = sqrt(diag(Pc));
sigmaX_RSS = computeRSS_coefErr(6, 26, 20, sigmaX);
sigmaC_RSS = computeRSS_coefErr(6, 26, 20, sigmaC);


%% COMPARISON WITH MC ANALYSIS

err_req = load("data/CoefErr_required.mat").CoefErr;
errRMS_req = load("data/CoefErr_RMS_requiered.mat").CoefErr_RMS;

% compute bounds. Px
boundUp = sigmaX' + median(err_req);
boundDown = -sigmaX' + median(err_req);

figure()
plot(1:46, err_req, 'o', 'Color','r', 'MarkerFaceColor','r')
hold on;
plot(1:46, boundUp, 1:46, boundDown, 'Color', 'g', 'LineWidth', 2)
title('coeficient error + Px bound')

% compute bounds. Px
boundUp = sigmaC' + median(err_req);
boundDown = -sigmaC' + median(err_req);

figure()
plot(1:46, err_req, 'o', 'Color','r', 'MarkerFaceColor','r')
hold on;
plot(1:46, boundUp, 1:46, boundDown, 'Color', 'b', 'LineWidth', 2)
title('coeficient error + Pc bound')


