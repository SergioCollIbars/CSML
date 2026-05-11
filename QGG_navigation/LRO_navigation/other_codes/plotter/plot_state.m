clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);
addpath('/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/other_codes/error_budget/');
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');

% colors
errorColor = [0.5 0.5 0.5]; % gray
boundColor = 'b';

%% PLOT POSITION & VELOCITY PERFORMANCE
% Author: Sergio Coll-Ibars
% Description: Plot state errors for the Low Lunar Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % folder = "/Users/sergiocollibars/Desktop/DLA_Gravity_Gradiometry/Simulator/data/";
folder = "/Users/sergiocollibars/OneDrive - UCB-O365/Kepler_codes/results/LRO_March_20_2015_12H/";

% set 4
% file = "orbit_LRO_ZZ_fc_1e-3_1mE_dataOut.mat";
file = "orbit_LRO_SIM2_dataOut.mat";

% load files
load(folder+file);

% time to UTC
utc  = cspice_et2utc(time', 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

% Monte Carlo runs
fields = fieldnames(mC_struct);
MC = numel(fields);
disp('Plot data ...');
figP = []; figV = [];
for mc  = 1:MC
    data = mC_struct.(fields{mc});

    % get rotation matrices
    BN_mat = compute_orientation_SC(time, data.Xf(1:6,:)', ...
        I_ALIG, NB_MOON_mat);

    % compute state error and covariance
    [errorP, errorV, sigmaP, sigmaV] = ...
            get_pos_vel_cov(data, length(time), BN_mat);
    
    % plot results
    if mc == 1
        figP = figure();
        for i = 1:3
            axP(i) = subplot(1,3,i);
            hold on
            grid on
            ylabel('Position Error')
        end
        figV = figure();
        for i = 1:3
            axV(i) = subplot(1,3,i);
            hold on
            grid on
            ylabel('Velocity Error')
        end
    end

    % Add covariance bounds
    if mc == MC
        figure(figP)
        for i = 1:3
            axes(axP(i))
            h1(i,:) = plot(tUTC,  3.*sigmaP(i,:), ...
                           tUTC, -3.*sigmaP(i,:), ...
                           'LineWidth', 2, ...
                           'Color', boundColor);
        end

        figure(figV)
        for i = 1:3
            axes(axV(i))
            h2(i,:) = plot(tUTC,  3.*sigmaV(i,:), ...
                           tUTC, -3.*sigmaV(i,:), ...
                           'LineWidth', 2, ...
                           'Color', boundColor);
        end
    end

    % Add Monte Carlo lines to position subplots
    figure(figP)
    for i = 1:3
        axes(axP(i))
        hP(i,mc) = plot(tUTC, errorP(i,:), ...
            'LineWidth', 1.0, 'Color', errorColor);
    end
    
    % Add Monte Carlo lines to velocity subplots
    figure(figV)
    for i = 1:3
        axes(axV(i))
        hV(i,mc) = plot(tUTC, errorV(i,:), ...
            'LineWidth', 1.0, 'Color', errorColor);
    end
end

% clear kernels
cspice_kclear

%% AUXILIARY FUNCTIONS
function [errorP, errorV, sigmaP, sigmaV] = ...
            get_pos_vel_cov(data, Nt, BN_mat)
    % state error: Inertial frame
    error  = data.state_true - data.Xf;
    
    sigmaP = nan(3, Nt); sigmaV = nan(3, Nt);
    errorP = nan(3, Nt); errorV = nan(3, Nt);
    for k = 1:Nt
        p = reshape(data.Pf(k, :), [18, 18]);
        maxInd = 3*k;
        minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);
    
        A   = blkdiag(BN, BN, eye(6),eye(6));
        P_B = A * p * A';
        std = sqrt(diag(P_B));
        sigmaP(:, k) = std(1:3);
        sigmaV(:, k) = std(4:6);
    
        errorP(:, k) = BN * error(1:3, k);
        errorV(:, k) = BN * error(4:6, k);
    end
end