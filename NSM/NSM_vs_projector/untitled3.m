clear
clc
close all

% Parameters
n_meas = 9;
n_x = 2;
n_c = 2;
n_trials = 30;
perturb_levels = logspace(0, 3, n_trials); % From small to large perturbations

rng(0); % Reproducibility

% True parameters
x_true = [1; -2];
c_true = [5; -6.0];

% Model matrices
A = randn(n_meas, n_x);
B = randn(n_meas, n_c);

% True noise covariance and weight
R_true = diag(linspace(0.5, 2.0, n_meas));
W_true = inv(R_true);

% Simulate true noise and observations
noise = mvnrnd(zeros(1, n_meas), R_true)';
l = A * x_true + B * c_true + noise;

% Left null space of B
N = null(B');

% Preallocate error arrays
err_nsm = zeros(1, n_trials);
err_proj = zeros(1, n_trials);

% Loop over perturbation levels
for i = 1:n_trials
    E = randn(n_meas);
    E = 0.5 * (E + E'); % Symmetrize

    % Perturbation matrix (symmetric noise)
    perturb = perturb_levels(i) * E ;
   % % perturb = diag([ones(1, n_meas)].*perturb_levels(i));

    % Perturbed weight matrix
    R_pert = R_true + perturb; 
    W_pert = inv(R_pert);

    % Ensure W_pert is positive definite (project to SPD if needed)
    [V,D] = eig(W_pert);
    D = max(D, 1e-6); % prevent negative eigenvalues
    W_pert = V * D * V';

    % NSM method
    A_nsm = N' * A;
    l_nsm = N' * l;
    W_nsm = N' * W_pert * N;
    x_hat_nsm = (A_nsm' * W_nsm * A_nsm) \ (A_nsm' * W_nsm * l_nsm);

    % Projector-based method
    P = eye(n_meas) - B * ((B' * W_pert * B) \ (B' * W_pert));
    A_proj = P * A;
    l_proj = P * l;
    x_hat_proj = (A_proj' * W_pert * A_proj) \ (A_proj' * W_pert * l_proj);

    % Store estimation errors
    err_nsm(i) = norm(x_hat_nsm - x_true);
    err_proj(i) = norm(x_hat_proj - x_true);
end

% Plotting
figure;
semilogx(perturb_levels, err_nsm, '-o', 'LineWidth', 2); hold on;
semilogx(perturb_levels, err_proj, '-s', 'LineWidth', 2);
xlabel('Perturbation Level ||\Delta W||');
ylabel('Estimation Error ||x_{est} - x_{true}||');
legend('NSM', 'Projector-Based');
title('Sensitivity to Weight Matrix Perturbation');
grid on;
