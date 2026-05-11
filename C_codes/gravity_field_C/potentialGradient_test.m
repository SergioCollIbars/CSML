clear;
clc; close all;
format long g;
addpath("MATLAB_functions/");

%% TEST MATLAB vs C++ potentialGradient_nm function
% Author: Sergio Coll-Ibars
% Date: 05/10/2026
% Description: Compare the potentialGradient_nm.m function to the C++ based
% routine. Evaluate computation speed and numerical accuracy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Planet constants (Moon)

GM = 4.902801076e12;         % [m^3/s^2]
R  = 1738e3;                 % [m]

r = (R + 50e3) * [1;1;1];   % [m]

[Cnm, Snm, ~] = readCoeff("HARMCOEFS_MOON_GRGM1200");

%% Speed and accuracy evaluation

n_array = 20:1200;

t_matlab = nan(size(n_array));
t_mex    = nan(size(n_array));

err_U    = nan(size(n_array));
err_grad = nan(size(n_array));
err_Hess = nan(size(n_array));

% Warm-up calls
potentialGradient_nm(Cnm, Snm, n_array(1), r, R, GM, 1);
potentialGradient_nm_mex(Cnm, Snm, n_array(1), r, R, GM, 1);

for idx = 1:length(n_array)

    n = n_array(idx);
    disp("Gravity degree: " + n);

    %% Check that both functions produce the same result
    [U_matlab, grad_matlab, Hess_matlab] = potentialGradient_nm( ...
        Cnm, Snm, n, r, R, GM, 1);

    [U_mex, grad_mex, Hess_mex] = potentialGradient_nm_mex( ...
        Cnm, Snm, n, r, R, GM, 1);

    err_U(idx) = abs(U_matlab - U_mex) / abs(U_matlab);

    err_grad(idx) = norm(grad_matlab - grad_mex) / ...
                    norm(grad_matlab);

    err_Hess(idx) = norm(Hess_matlab - Hess_mex, 'fro') / ...
                    norm(Hess_matlab, 'fro');

    %% Time both functions
    f_matlab = @() potentialGradient_nm(Cnm, Snm, n, r, R, GM, 1);
    f_mex    = @() potentialGradient_nm_mex(Cnm, Snm, n, r, R, GM, 1);

    t_matlab(idx) = timeit(f_matlab);
    t_mex(idx)    = timeit(f_mex);

end

%% Display maximum relative differences

fprintf('\nMaximum relative differences:\n');
fprintf('Potential:        %.3e\n', max(err_U));
fprintf('First gradient:   %.3e\n', max(err_grad));
fprintf('Second gradient:  %.3e\n', max(err_Hess));

%% Plot absolute execution times

figure();
plot(n_array, t_matlab, 'LineWidth', 2); hold on;
plot(n_array, t_mex, 'LineWidth', 2);
grid on;
xlabel('Maximum degree');
ylabel('Execution time [s]');
legend('MATLAB', 'MEX', 'Location', 'northwest');

%% Plot speedup

figure();
plot(n_array, t_matlab ./ t_mex, 'LineWidth', 2);
grid on;
xlabel('Maximum degree');
ylabel('Speedup: MATLAB / MEX');

%% Plot relative differences

figure();
plot(n_array, err_U, 'LineWidth', 2); hold on;
plot(n_array, err_grad, 'LineWidth', 2);
plot(n_array, err_Hess, 'LineWidth', 2);
grid on;
xlabel('Maximum degree');
ylabel('Relative difference');
legend('Potential', 'First gradient', 'Second gradient', ...
    'Location', 'best');