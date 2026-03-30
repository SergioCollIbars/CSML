%% demo_projection_vs_coestimation_NOnoise.m
% Deterministic comparison:
%  (A) Co-estimation of [x; eta]
%  (B) Projection with P_G = I - G(G'G)^{-1}G'
%
% No noise, no weighting, no Schur complements.

clear; clc; close all;
rng(4);

%% Dimensions
m  = 120;    % number of measurements
nx = 6;      % target parameters
ng = 4;      % nuisance parameters

%% Construct deterministic matrices
H = randn(m,nx);
G = randn(m,ng);

%% Choose arbitrary but consistent truth
x_true   = randn(nx,1);
eta_true = randn(ng,1);

% Exact measurements (no noise)
y = H*x_true + G*eta_true;

%% ===============================
% (A) Co-estimation
%% ===============================
A = [H G];

% Least-squares solution
theta_hat = (A' * A) \ (A' * y);
x_hat_co  = theta_hat(1:nx);

% Formal covariance (up to scale factor)
Cov_full = inv(A' * A);
Cov_x_co = Cov_full(1:nx,1:nx);

%% ===============================
% (B) Projection approach
%% ===============================
PG = eye(m) - G * ((G' * G) \ G');

% Projected LS
x_hat_proj = (H' * PG * H) \ (H' * PG * y);

% Formal covariance
Cov_x_proj = inv(H' * PG * H);

%% ===============================
% Comparisons
%% ===============================
fprintf('||x_hat_co - x_hat_proj||_2                  = %.3e\n', ...
        norm(x_hat_co - x_hat_proj));

fprintf('||Cov_x_co - Cov_x_proj||_F                  = %.3e\n', ...
        norm(Cov_x_co - Cov_x_proj,'fro'));

fprintf('Relative x error                             = %.3e\n', ...
        norm(x_hat_co - x_hat_proj)/max(1,norm(x_hat_co)));

fprintf('Relative covariance error                   = %.3e\n', ...
        norm(Cov_x_co - Cov_x_proj,'fro')/max(1,norm(Cov_x_co,'fro')));

%% Optional: verify projector properties
fprintf('||PG - PG''||_F (symmetry)                    = %.3e\n', ...
        norm(PG - PG','fro'));

fprintf('||PG^2 - PG||_F (idempotent)                 = %.3e\n', ...
        norm(PG*PG - PG,'fro'));
