clear;
clc;
close all;
format long g;

%%                         NAVIGATION TEST                               %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/31/2024                                                      %
%                                                                         %
%   Description: Code to test gradiometer ability to reconstruct          %
%   trajectory based on the QGG measurements. Read trajectory from QGG    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultAxesFontSize',16);
addpath("functions/")
addpath('../functions/');
addpath('../dataFiles/');
addpath('../../../QGG/data_files/');

dOrbit = readtable("orbitData.txt");
dAcc = readtable('accData.txt');
dataOrbit = table2array(readtable("orbitData.txt"));
dataMeas = table2array(readtable('accData.txt'));

% number of points
t = dataOrbit(:, 1);

% Planet conditions
GM = 5.2;       % [m^3/s^2]
n_max = 9;
normalized = 1;
path = "HARMCOEFS_BENNU_CD_1.txt";
[Cmat, Smat, Re] = readCoeff(path);
planetParams = [GM, Re, n_max, normalized];

W = 4.06130329511851E-4;          % [rad/s]
W0 = 0;                           % [rad/s]
RA = deg2rad(86.6388);            % [rad]
DEC = deg2rad(-65.1086);          % [rad]
poleParams = [W, W0, RA, DEC];

% truth position
rt_ACI = dataOrbit(:, 2:4)';
vt_ACI = [dataOrbit(:, 8), dataOrbit(:, 9), dataOrbit(:, 10)]';
state_true = [rt_ACI; vt_ACI]'; 

% Nominal orbit. Initial conditions
X0 = [rt_ACI(:, 1); vt_ACI(:, 1)];
delta = [rt_ACI(:, 2); vt_ACI(:, 2)] -X0;

% generate gradiometer measurements
[ddU_ACI] = gradiometer_meas(t, state_true(:, 1:3), Cmat, Smat, ...
    poleParams, planetParams);

% reconstruct trajectory
[delta_new] = reconstruct_traj(t, ddU_ACI, delta);

 % state reconstruction from nominal
state_recons = zeros(6, length(t));
state_recons(:, 1) = X0;
for j = 2:length(t)
    state_recons(:, j) = state_recons(:, j-1) + delta_new(:, j);
end

%% PLOTS

% Plot trajectory in Inertial frame
figure();
plot3(state_true(:, 1)', state_true(:, 2)', state_true(:, 3)', ...
    'LineWidth', 2, 'Color', 'b');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Real trajectory vs Nominal trajectory. Inertial Frame');
axis equal;
grid on;

left = [1, 3, 5];
right = [2, 4, 6];

% % % plot true and nominal trajectory
% % figure()
% % for j =1:3
% %     subplot(3, 2, left(j));
% %     plot(t./86400, state_true(:, j), 'LineWidth', 2, 'Color', 'b')
% %     hold on;
% %     plot(t./86400, state_nom(:, j), 'LineWidth', 2, 'Color', 'r')
% %     xlabel('TIME [days]')
% %     ylabel('r_' + string(j))
% %     legend('true', 'nominal')
% % 
% %     subplot(3, 2, right(j));
% %     plot(t./86400, state_true(:, j) - state_nom(:, j), ...
% %         'LineWidth', 2, 'Color', "#D95319")
% %     xlabel('TIME [days]')
% %     ylabel('error [m]')
% % end
% % sgtitle('True and nominal trajectory  components. Inertial Frame')

% plot reconstructed trajectory
figure()
for j =1:3
    subplot(3, 2, left(j));
    plot(t./86400, state_true(:, j), 'LineWidth', 2, 'Color', 'b')
    hold on;
    plot(t./86400, state_recons(j, :), 'LineWidth', 2, 'Color', 'g')
    xlabel('TIME [days]')
    ylabel('r_' + string(j))
    legend('true', 'reconstructed')

    subplot(3, 2, right(j));
    plot(t./86400, state_true(:, j) - state_recons(j, :)', ...
        'LineWidth', 2, 'Color', "#FF00FF")
    xlabel('TIME [days]')
    ylabel('error [m]')
end
sgtitle('Reconstructed trajectory components. Inertial Frame')

% % % plot true and nominal velocity
% % figure()
% % for j =1:3
% %     subplot(3, 2, left(j));
% %     plot(t./86400, state_true(:, j+3), 'LineWidth', 2, 'Color', 'b')
% %     hold on;
% %     plot(t./86400, state_nom(:, j+3), 'LineWidth', 2, 'Color', 'r')
% %     xlabel('TIME [days]')
% %     ylabel('v_' + string(j))
% %     legend('true', 'nominal')
% % 
% %     subplot(3, 2, right(j));
% %     plot(t./86400, state_true(:, j+3) - state_nom(:, j+3), ...
% %         'LineWidth', 2, 'Color', "#D95319")
% %     xlabel('TIME [days]')
% %     ylabel('error [m]')
% % end
% % sgtitle('True and nominal velocity  components. Inertial Frame')

% plot reconstructed velocity
figure()
for j =1:3
    subplot(3, 2, left(j));
    plot(t./86400, state_true(:, j+3), 'LineWidth', 2, 'Color', 'b')
    hold on;
    plot(t./86400, state_recons(j+3, :), 'LineWidth', 2, 'Color', 'g')
    xlabel('TIME [days]')
    ylabel('v_' + string(j))
    legend('true', 'reconstructed')

    subplot(3, 2, right(j));
    plot(t./86400, state_true(:, j+3) - state_recons(j+3, :)', ...
        'LineWidth', 2, 'Color', "#FF00FF")
    xlabel('TIME [days]')
    ylabel('error [m]')
end
sgtitle('Reconstructed velocity components. Inertial Frame')

% % 
% % % plot showing error comparison velocity
% % figure()
% % for j =1:3
% %     subplot(3, 2, left(j));
% %     plot(t./86400, abs(state_true(:, j) - state_nom(:, j)),...
% %         'LineWidth', 2, 'Color', "#D95319")
% %     hold on;
% %     plot(t./86400, abs(state_true(:, j) - state_recons(j, :)'),...
% %         'LineWidth', 2, 'Color', "#FF00FF")
% %     xlabel('TIME [days]')
% %     ylabel('r_' + string(j))
% %     legend('true - nom', 'true - recons')
% % 
% %     subplot(3, 2, right(j));
% %     plot(t./86400, abs(state_true(:, j+3) - state_nom(:, j+3)), ...
% %         'LineWidth', 2, 'Color', "#D95319")
% %     hold on;
% %     plot(t./86400, abs(state_true(:, j+3) - state_recons(j+3, :)'), ...
% %         'LineWidth', 2, 'Color', "#FF00FF")
% %    xlabel('TIME [days]')
% %    ylabel('v_' + string(j))
% %    legend('true - nom', 'true - recons')
% % end
% % sgtitle('Absolute trajectory and velocity error components. Inertial Frame')

% % % plot state deviations
% % figure()
% % for j =1:3
% %     subplot(3, 2, left(j));
% %     plot(t./86400, delta_new(j, :),...
% %         'LineWidth', 2, 'Color', "k")
% %     xlabel('TIME [days]')
% %     ylabel('r_' + string(j))
% %   
% % 
% %     subplot(3, 2, right(j));
% %     plot(t./86400, delta_new(j+3, :), ...
% %         'LineWidth', 2, 'Color', "k")
% %    xlabel('TIME [days]')
% %    ylabel('v_' + string(j))
% % end
% % sgtitle('trajectory and velocity state deviations. Inertial Frame')


