clc;
close all;
clear;
format long g;

%%             POTENTIAL RECURSION FORMULATION TEST

% Decription: Given an orbit, compute the potential gradient (acc) using
% the recursive formulation.

% Imports
addpath("functions/");
addpath("//Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/")

% Inputs
N = 20;                             % Number of points

phi = linspace(0, 2*pi, N);         % Latitude [rad]
%lambda = linspace(0, 2*pi, N);      % Longitude [rad]
lambda = ones(1, N) * 0;

rn = 1000;                          % orbit radius [m]
GM = 5.2;                           % gravity param [m^3 / s^2]

% Hamonics coefficients
path = "HARMCOEFS_BENNU_CD_1.txt";
[C, S, Re, normalized] = readCoeff(path);
path = "HARMCOEFS_BENNU_CD_0.txt";
[C2, S2, ~, ~] = readCoeff(path);

n_max = length(C) - 1;            % spherical harmonic order

% Orbit position
x = rn*cos(phi).*cos(lambda);
y = rn*cos(phi).*sin(lambda);
z = rn*sin(phi);

% Orbit radius vector. Inertial frame
r = [x; y; z];
v = [0; sqrt(GM/rn); 0];

% gravity partials vector definition
U_err = zeros(1, N);
dU_err = zeros(3, N);
ddU_err = zeros(3*N, 3);
acc = zeros(3, N);

% error vector and tolerance
err_relative = ones(9, N);
tolerance = 1E-6;

% spherical harmonics function options
option = ["numerical", "cartesian", "0"];

% Numerical vs analytical computation
disp('Starting test...')
for k = 1:N
    disp("  iter: " + string(k));
    % rotate form RTN to ECI
    [NB] = RTN2ECI(r(:, k), v);
    BN = NB';

    % current position
    rk = r(:, k);                               % ECI coords
    rk_RTN = BN * rk;                           % ECEF coords
    phik = phi(k);
    lambdak = lambda(k);
    rnk = vecnorm(rk);

    % compute spherical partials. Visibility matrix
    [H] = potentialGradient_Cnm(n_max, rk_RTN, Re, GM, eye(3), normalized);

    % compute gravity potential gradient. Numerical
    [U_n, dU_n, ddU_n] = potentialGradient_nm(C, S, n_max, rk_RTN, ...
        Re, GM, normalized);

    % compute gravity potential. Analytically
    [U_s, dU_s, ddU_s] = sphericalHarmonics(GM, Re, n_max, C2, S2, ...
        rk_RTN, option);

    % save truth value acc
    acc(:, k) = dU_s;
    
    % compute error
    analytical = reshape(ddU_s, [9, 1]);
    numerical = reshape(ddU_n, [9, 1]);

    for j = 1:9
        if(analytical(j) == 0)
            err_relative(j, k) = abs(analytical(j) - numerical(j));
        else
            err_relative(j, k) = abs((analytical(j) - numerical(j)) /...
                analytical(j));
        end
    end
end

% print test results
maxVal = max(err_relative');
disp('Relative error per component: ' + string(maxVal'))
if(maxVal < tolerance)
    disp('  TEST PASS SUCCESFULLY');
else
    disp('  TEST FAILED');
end






