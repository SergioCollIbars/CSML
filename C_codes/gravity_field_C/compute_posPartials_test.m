clear;
clc;
close all;
format long g;

addpath("MATLAB_functions/");
addpath(genpath("functions/"));

%% TEST MATLAB vs C++ compute_posPartials function
% Author: Sergio Coll-Ibars
% Date: 05/11/2026
%
% Description:
% Compare the analytical MATLAB compute_posPartials.m function against the
% analytical MEX implementation compute_posPartials_analytic_mex.cpp.
%
% The script checks:
%   1. numerical agreement
%   2. execution time
%   3. speedup factor

%% Planet constants: Moon

GM = 4.902801076e12;     % [m^3/s^2]
R  = 1738e3;             % [m]

%% Test position

r = (R + 50e3) * [1; 1; 1];   % [m]

%% Attitude / rotation matrices

% For a first test, use identity rotations.
% ACAF_ACI maps inertial coordinates to body-fixed coordinates.
% ACAF_B maps body-frame gradiometer components to the body-fixed frame,
% depending on your convention inside compute_posPartials.
ACAF_ACI = eye(3);
ACAF_B   = eye(3);

%% Load gravity coefficients

[Cnm, Snm, ~] = readCoeff("HARMCOEFS_MOON_GRGM1200");

%% Test settings

normalized = 1;

n_array = 20:1:1200;

t_matlab = nan(size(n_array));
t_mex    = nan(size(n_array));

err_abs = nan(size(n_array));
err_rel = nan(size(n_array));

%% Warm-up calls

n0 = n_array(1);

H_matlab = compute_posPartials( ...
    n0, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

H_mex = compute_posPartials_mex( ...
    n0, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

fprintf("Warm-up relative error: %.3e\n", ...
    norm(H_matlab - H_mex, "fro") / norm(H_matlab, "fro"));

%% Speed and accuracy evaluation

for idx = 1:numel(n_array)

    n = n_array(idx);

    fprintf("Gravity degree: %d\n", n);

    %% Accuracy check

    H_matlab = compute_posPartials( ...
        n, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

    H_mex = compute_posPartials_mex( ...
        n, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

    err_abs(idx) = max(abs(H_matlab(:) - H_mex(:)));

    err_rel(idx) = norm(H_matlab - H_mex, "fro") / ...
                   max(norm(H_matlab, "fro"), eps);

    %% Timing

    f_matlab = @() compute_posPartials( ...
        n, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

    f_mex = @() compute_posPartials_mex( ...
        n, normalized, Cnm, Snm, R, GM, r, ACAF_ACI, ACAF_B);

    t_matlab(idx) = timeit(f_matlab);
    t_mex(idx)    = timeit(f_mex);

end

%% Display maximum differences

fprintf("\nMaximum differences:\n");
fprintf("Max absolute error: %.3e\n", max(err_abs));
fprintf("Max relative error: %.3e\n", max(err_rel));

fprintf("\nTiming summary:\n");
fprintf("Fastest MATLAB time: %.3e s\n", min(t_matlab));
fprintf("Fastest MEX time:    %.3e s\n", min(t_mex));
fprintf("Max speedup:         %.2f x\n", max(t_matlab ./ t_mex));

%% Plot absolute execution times

figure();
plot(n_array, t_matlab, "LineWidth", 2);
hold on;
plot(n_array, t_mex, "LineWidth", 2);
grid on;
xlabel("Maximum degree");
ylabel("Execution time [s]");
legend("MATLAB", "MEX", "Location", "northwest");
title("compute\\_posPartials Execution Time");

%% Plot speedup

figure();
plot(n_array, t_matlab ./ t_mex, "LineWidth", 2);
grid on;
xlabel("Maximum degree");
ylabel("Speedup: MATLAB / MEX");
title("compute\\_posPartials MEX Speedup");

%% Plot relative differences

figure();
semilogy(n_array, err_rel, "LineWidth", 2);
grid on;
xlabel("Maximum degree");
ylabel("Relative difference");
title("MATLAB vs MEX Relative Difference");

%% Plot absolute differences

figure();
semilogy(n_array, err_abs, "LineWidth", 2);
grid on;
xlabel("Maximum degree");
ylabel("Maximum absolute difference");
title("MATLAB vs MEX Absolute Difference");