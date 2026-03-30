clear;
clc;
close all;


% Parameters
rng(0);  % For reproducibility
m = 5;         % Original dimension
r = 3;         % Rank of A (number of constraints)
nSamples = 100000;

% Step 1: Diagonal covariance matrix R (different variances)
std_vec = [1.0, 2.0, 1.5, 0.8, 1.2];  % Standard deviations
% % std_vec = [1, 1, 1, 1, 1]; 
R = diag(std_vec .^ 2);              % Covariance matrix (diagonal)

% Generate Gaussian samples with covariance R
L = chol(R, 'lower');                % Cholesky decomposition
samples = L * randn(m, nSamples);    % Each column is a sample

% Step 2: Define matrix A and compute its left null space N
A = randn(m, r);                     % Constraint matrix
N = null(A');                   % Left null space basis (m x (m - r))
V = N * N';

f1 = A * inv(A' * (R\A)) *(A'/(R));
F = eye(m,m) - f1;

% Step 3: Project samples
projected_samples = N' * samples;    % Size: (m - r) x nSamples
projected_samples = V * samples;    % Size: (m - r) x nSamples

% Step 4: Sample covariance of projected data
sample_cov_proj = cov(projected_samples');

% Step 5: Theoretical matrices
theoretical1 = N' * inv(R) * N;
theoretical1 = V * inv(R) * V';
theoretical2 = inv(N' * R * N);

% Step 6: Compute norms of differences
diff1 = norm(inv(sample_cov_proj) - theoretical1);
diff2 = norm(inv(sample_cov_proj) - theoretical2);

% Display results
fprintf('Difference from N^T * inv(R) * N: %.6f\n', diff1);
fprintf('Difference from inv(N^T * R * N): %.6f\n', diff2);

% Optional: visualize the projected samples (only if m-r = 2)
if size(projected_samples,1) == 2
    figure;
    scatter(projected_samples(1,:), projected_samples(2,:), 10, 'filled');
    axis equal;
    title('Projected Gaussian Cloud (2D)');
end
