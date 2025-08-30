clear;
clc;
% Dimensions
m = 6;     % total observations
n = 3;     % dimension of H
p = 2;     % dimension of B

epsilon = 1e-6;

% Random full-rank B and H matrices
B = randn(m, p);
H = randn(m, n);

% Compute null space of B'
C = null(B');

% Ensure R is diagonal and positive
R_diag = 1 + rand(m, 1);

% Compute projected covariance and w = inv(C' * R * C)
A = C' * diag(R_diag) * C;
w = inv(A);

% Compute main expression
M = (C'*H)' * w * (C'*H);

% Allocate outputs
dM_analytic = cell(m, 1);
dM_numeric = cell(m, 1);

% Loop over each diagonal element r_i
for i = 1:m
    % ---- Analytical ----
    ci = C(i, :)';                   % (m - rank(B)) x 1
    dA = ci * ci';                   % outer product
    dM_analytic{i} = -(C'*H)' * w * dA * w * (C'*H);

    % ---- Numerical ----
    R_perturbed = R_diag;
    R_perturbed(i) = R_perturbed(i) + epsilon;
    A_eps = C' * diag(R_perturbed) * C;
    w_eps = inv(A_eps);
    M_eps = (C'*H)' * w_eps * (C'* H);

    dM_numeric{i} = (M_eps - M) / epsilon;

    % ---- Comparison ----
    fprintf('i = %d | norm diff: %.2e | rel diff: %.2e\n', ...
        i, norm(dM_numeric{i} - dM_analytic{i}), ...
        norm(dM_numeric{i} - dM_analytic{i}) / norm(dM_analytic{i}));
end
