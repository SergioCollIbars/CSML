clear;
clc;
close all;

% Dimensions
n = 9;   % number of measurements
p = 3;   % number of estimated parameters
q = 3;   % number of constraints (so null space dim = n - q)

% Synthetic data
H = randn(n, p);
y = randn(n, 1);
B = randn(n, q);
C = null(B');  % C is n x (n - q)

% Original covariance
R = diag(1 + 0.1*rand(n,1));  % nominal diagonal covariance
R_proj = C' * R * C;

% Projected matrices
Hhat = C' * H;
yhat = C' * y;

% LS solution
A = Hhat' * (R_proj \ Hhat);
b = Hhat' * (R_proj \ yhat);
x_nom = A \ b;

% Perturb R
epsilon = 1e-6;
deltaR = diag(randn(n,1));  % perturbation (diagonal for simplicity)

% Projected perturbation
deltaR_proj = C' * deltaR * C;

% Analytical sensitivity
dA = -Hhat' * (R_proj \ deltaR_proj / R_proj) * Hhat;
db = -Hhat' * (R_proj \ deltaR_proj / R_proj) * yhat;
dx_analytic = -A \ dA * x_nom + A \ db;

% Finite difference
R_pert = R + epsilon * deltaR;
R_proj_pert = C' * R_pert * C;
A_pert = Hhat' * (R_proj_pert \ Hhat);
b_pert = Hhat' * (R_proj_pert \ yhat);
x_pert = A_pert \ b_pert;

dx_fd = (x_pert - x_nom) / epsilon;

% Display error
disp('Analytical dx:');
disp(dx_analytic);

disp('Finite-diff dx:');
disp(dx_fd);

disp('Error norm:');
disp(norm(dx_analytic - dx_fd));

%%

% Nominal covariance
invR = inv(R);

% Define projector P
BtRB_inv = inv(B' * invR * B);
P = eye(n) - B * BtRB_inv * B' * invR;

% Projected matrices
Hhat = P * H;
yhat = P * y;

% LS solution
A = Hhat' * invR * Hhat;
b = Hhat' * invR * yhat;
x_nom = A \ b;

% Perturbed covariance
R_pert = R + epsilon * deltaR;
invR_pert = inv(R_pert);

% Perturbed projector
BtRB_pert_inv = inv(B' * invR_pert * B);
P_pert = eye(n) - B * BtRB_pert_inv * B' * invR_pert;

% Perturbed projection
Hhat_pert = P_pert * H;
yhat_pert = P_pert * y;

% Perturbed LS solution
A_pert = Hhat_pert' * invR_pert * Hhat_pert;
b_pert = Hhat_pert' * invR_pert * yhat_pert;
x_pert = A_pert \ b_pert;

% Finite difference sensitivity
dx_fd = (x_pert - x_nom) / epsilon;

disp('Finite-difference dx:');
disp(dx_fd);

