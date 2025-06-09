clear;
clc;
close all;

format long g;
addpath('../NSM/functions/')
addpath("config/")
addpath("data_files/")
set(0,'defaultAxesFontSize',16);
%%              ATTITUDE SENSITIVITY ANALYSIS
% Date: 05/27/2025
% Author: Sergio Coll Ibars
% Description: given a trajectory, compute the sensitivity of gradiometer
% measurements to a fix attitude error.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load trajectory data
T = readtable("orbitData.txt");
t  = T.TIME;
r_ACI = [T.ri_x, T.ri_y, T.ri_z]';   % [m]
v_ACI = [T.vi_x, T.vi_y, T.vi_z]';   % [m/s]
phase = "B";    % options: A or B

% get BENNU coefficients
path = "HARMCOEFS_BENNU_OSIRIS_0.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.19;
n_max  = 6;
normalized = 0;
asterParams = [GM, Re, n_max, normalized];
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

[Nc, Ns, Ncs] = count_num_coeff(n_max); 
X       = mat2list(Cnm, Snm, Nc, Ns);

% noise parameters
noise0 = zeros(9, 1);
sigma_ii = 6.32E-12 * 10;
sigma_ij = 2.52E-10 * 10;
R        = diag([sigma_ii, sigma_ij, sigma_ij, sigma_ii, sigma_ij, sigma_ii].^2); 

% attitude errors
dRA  = normrnd(0, deg2rad(0.007),  [1, length(t)]);
dDEC = normrnd(0, deg2rad(0.003),  [1, length(t)]);
dW   = normrnd(0, deg2rad(0.012)/86400, [1, length(t)]);

% Compute sensitivity error
std_theta = 10 * pi/(180*3600);                          % [rad]
% % A_theta   = normrnd(0, std_theta,  [3, length(t)]);      % [rad]
A_theta   = std_theta.*ones(3, length(t));
S         = ones(9, length(t)) * NaN; 
Err       = S; Err_SC = S;
Ax = 0; S_planet = 0; S_SC = 0;
for j = 1:length(t)
    disp(j/length(t));
    
    % planet orientation
    Wt  = W * t(j) + W0;
    dWT = (W + dW(j)) * t(j) + W0;
    ACAF_ACI_true = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    ACAF_ACI_nom  = rotationMatrix(pi/2 + RA + dRA(j), pi/2 - (DEC + dDEC(j)), dWT, [3, 1, 3]);

    % compute gradiometer partials
    ACI_BODY = eye(3,3);
    [Y_true, ~, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI_true, [r_ACI(:, j)', v_ACI(:, j)'], ...
            noise0, Cnm, Snm, ACI_BODY);
    [Y_planet, H_planet_ACI, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI_nom, [r_ACI(:, j)', v_ACI(:, j)'], ...
            noise0, Cnm, Snm, ACI_BODY);
    
    ACI_BODY = rotationMatrix(A_theta(1, j), A_theta(2, j), A_theta(3, j), [3, 2, 1]);
    [~, ~, H_SC_BODY] = gradiometer_meas(t(j) ,[GM, Re, n_max, normalized], ACAF_ACI_true, [r_ACI(:, j)', v_ACI(:, j)'], ...
            noise0, Cnm, Snm, ACI_BODY);
    
    % rotation to SC body frame
    ACI_BODY = rotationMatrix(A_theta(1, j), A_theta(2, j), A_theta(3, j), [3, 2, 1]);
    T_true = [Y_true(1), Y_true(2), Y_true(3);...
        Y_true(4), Y_true(5), Y_true(6);...
        Y_true(7), Y_true(8), Y_true(9)];
    T_SC = ACI_BODY' * (T_true) * ACI_BODY;
    Y_SC = [T_SC(1,1); T_SC(1,2);T_SC(1,3);T_SC(2,1);T_SC(2,2);T_SC(2,3);...
        T_SC(3,1);T_SC(3,2);T_SC(3,3)];

    Err_SC(:, j) = (Y_true) - Y_SC;

    % compute rotation partials
    [Hrot]  = compute_rotPartials_analy(Y_true);
    S(:, j) = Hrot * A_theta(:, j);

    % compute error due to planet rotation
    Err(:, j) = Y_true - Y_planet;

    % visibility matrix
    H = [H_planet_ACI(1, :); H_planet_ACI(4, :); H_planet_ACI(7, :);H_planet_ACI(5, :);...
        H_planet_ACI(8, :); H_planet_ACI(9, :)];

    % planet orientation sensitivity
    S_planet = S_planet + H' * inv(R) * [Err(1, j);Err(2, j);Err(3,j);Err(5,j);Err(6,j);Err(9, j)];

    % S/C orientation sensitivity
    H = [H_SC_BODY(1, :); H_SC_BODY(4, :); H_SC_BODY(7, :);H_SC_BODY(5, :);...
        H_SC_BODY(8, :); H_SC_BODY(9, :)];

    S_SC = S_SC +  H' * inv(R) * [Err_SC(1, j);Err_SC(2, j);Err_SC(3,j);Err_SC(5,j);Err_SC(6,j);Err_SC(9, j)];

    % information matrix
    Ax = Ax + H' * inv(R) * H;
end

% STD
sigma = sqrt(diag(inv(Ax)));
SH_err_SC = inv(Ax) * S_SC;

SH_err_planet = inv(Ax) * S_planet;

sigma_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                sigma, zeros(n_max+1, n_max+1), ...
                zeros(n_max+1, n_max+1)); 

RMS_truth  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                X, zeros(n_max+1, n_max+1), ...
                zeros(n_max+1, n_max+1)); 

% SH error due to planet orientation
SH_err_planet  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                SH_err_planet, zeros(n_max+1, n_max+1), ...
                zeros(n_max+1, n_max+1)); 

% SH error due to SC orientation
SH_err_SC  = computeRMS_coeffErr(n_max, Nc, Ns, ...
                SH_err_SC, zeros(n_max+1, n_max+1), ...
                zeros(n_max+1, n_max+1)); 

figure()
semilogy(1:n_max, RMS_truth, 'LineWidth', 2, 'Color','k')
hold all;
semilogy(1:n_max, sigma_RMS, 'LineWidth', 2, 'Color','b')
semilogy(1:n_max, SH_err_planet, 'LineWidth', 2, 'Color','r')
xlabel('SH degree [n]')
title('SH RMS error due to planet orientation errors')
legend('thuth', '1-\sigma', 'RMS error')

figure()
semilogy(1:n_max, RMS_truth, 'LineWidth', 2, 'Color','k')
hold all;
semilogy(1:n_max, sigma_RMS, 'LineWidth', 2, 'Color','b')
semilogy(1:n_max, SH_err_SC, 'LineWidth', 2, 'Color','r')
xlabel('SH degree [n]')
title('SH RMS error due to S/C orientation errors')
legend('thuth', '1-\sigma', 'RMS error')


% plot measurement residuals
figure()
plot(t./86400, [S(1, :);S(2, :);S(3,:);S(5,:);S(6,:);S(9, :)]./1E-12, 'LineWidth', 2)
legend('\Gamma_{xx}', '\Gamma_{xy}', '\Gamma_{xz}', '\Gamma_{yy}', ...
    '\Gamma_{yz}', '\Gamma_{zz}');
title('Measurement residuals due to S/C attitude errors')
xlabel('TIME [days]')
ylabel('milli-Eotvos [mE]')

figure()
plot(t./86400, [Err(1, :);Err(2, :);Err(3,:);Err(5,:);Err(6,:);Err(9, :)]./1E-12, 'LineWidth', 2)
legend('\Gamma_{xx}', '\Gamma_{xy}', '\Gamma_{xz}', '\Gamma_{yy}', ...
    '\Gamma_{yz}', '\Gamma_{zz}')
title('Measurement residuals due to asteroid attitude errors')
xlabel('TIME [days]')
ylabel('milli-Eotvos [mE]')

%% FUNCTIONS
function [Hrot] = compute_rotPartials_analy(Y)
    Yxx = Y(1); Yxy = Y(2); Yxz = Y(3);
    Yyy = Y(5); Yyz = Y(6); Yzz = Y(9);

    Hrot = -[0, 2*Yxz, -2*Yxy;...
           -Yxz, Yyz, Yxx - Yyy;...
        Yxy, Yzz - Yxx, -Yyz;...
        -Yxz, Yyz, Yxx - Yyy;...
        -2*Yyz, 0, 2*Yxy;...
        Yyy - Yzz, -Yxy, Yxz;...
        Yxy, Yzz - Yxx, -Yyz;...
        Yyy - Yzz, -Yxy, Yxz;...
        2*Yyz, -2*Yxz, 0];
end