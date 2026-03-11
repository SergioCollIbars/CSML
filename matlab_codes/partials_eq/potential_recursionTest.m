clc;
close all;
clear;
format long g;

%%             POTENTIAL RECURSION FORMULATION TEST

% Decription: Given an orbit, compute the potential gradient (acc) using
% the recursive formulation.

% Imports
addpath("functions/");

% Inputs
N = 20;                             % Number of points

phi = linspace(0, 2*pi, N);         % Latitude [rad]
lambda = linspace(0, 2*pi, N);      % Longitude [rad]
%lambda = ones(1, N) * 0;

rn = 500;                          % orbit radius [m]
GM = 5.2;                          % gravity param [m^3 / s^2]


% Hamonics coefficients
path = "/Users/sergiocollibars/Desktop/CSML/QGG_gravEstim/data_files/HARMCOEFS_BENNU_CD_1.txt";
[C, S, Re, normalized] = readCoeff(path);
path = "/Users/sergiocollibars/Desktop/CSML/QGG_gravEstim/data_files/HARMCOEFS_BENNU_CD_0.txt";
[C2, S2, ~, ~] = readCoeff(path);

n_max = length(C) - 1;            % spherical harmonic order
n_max = 0;

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

C3 = C.*0;
S3 = S.*0;
C3(1,1)  = 1;

% spherical harmonics function options
option = ["numerical", "cartesian", "1"];

% Numerical vs analytical computation
disp('Starting test...')
for k = 1:N
    disp("  iter: " + string(k));
    % rotate form RTN to ECI
    BN = eye(3);

    % current position
    rk = r(:, k);                               % ECI coords
    rk_RTN = BN * rk;                           % ECEF coords
    phik = phi(k);
    lambdak = lambda(k);
    rnk = vecnorm(rk);

    % compute spherical partials. Visibility matrix
    [H] = potentialGradient_Cnm(n_max, rk_RTN, Re, GM, eye(3), normalized);

    % compute gravity potential gradient. Numerical
    [U_n, dU_n, ddU_n] = potentialGradient_nm(C3, S3, n_max, rk_RTN, ...
        Re, GM, normalized);
    
    % version 2
    % [U, a, T] = gravity_potential_acc_tensor(C, S, rk_RTN, Re, GM, normalized);

    % compute gravity potential. Analytically
%     [U_s, dU_s, ddU_s] = sphericalHarmonics(GM, Re, n_max, C, S, ...
%         rk_RTN, option);

    [U_val, a_val, ddU_s, ~, ~, ~] = ...
    gravity_symbolic_eval(C3, S3, rk_RTN, Re, GM, normalized);

    % save truth value acc
    acc(:, k) = a_val;
    
    % compute error
    analytical = reshape(ddU_s, [9, 1]);
    numerical = reshape(ddU_n, [9, 1]);

    for j = 1:9
        if(analytical(j) == 0)
            err_relative(j, k) = abs(analytical(j) - numerical(j));
        else
            err_relative(j, k) = abs((analytical(j) - numerical(j))./analytical(j));
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


%% FUNCTIONS

function [U_val, a_val, T_val, U_sym, a_sym, T_sym] = ...
    gravity_symbolic_eval(Cnm, Snm, r_ecef, Re, GM, normalized)
%GRAVITY_SYMBOLIC_EVAL
% Build the gravity potential symbolically and evaluate:
%   - potential U
%   - acceleration a = grad(U)
%   - gravity tensor T = Hessian(U)
%
% Inputs
%   Cnm        : (N+1)x(N+1) cosine coefficient matrix, Cnm(n+1,m+1)
%   Snm        : (N+1)x(N+1) sine   coefficient matrix, Snm(n+1,m+1)
%   r_ecef     : 3x1 Cartesian evaluation point [x;y;z]
%   Re         : reference radius
%   GM         : gravitational parameter
%   normalized : true for fully-normalized coefficients, false otherwise
%
% Outputs
%   U_val      : potential evaluated at r_ecef
%   a_val      : 3x1 acceleration evaluated at r_ecef
%   T_val      : 3x3 gravity tensor evaluated at r_ecef
%   U_sym      : symbolic potential U(x,y,z)
%   a_sym      : symbolic gradient
%   T_sym      : symbolic Hessian
%
% Notes
%   - This is symbolic, so it gets slow quickly as degree increases.
%     It is practical only for modest n_max.
%   - Coefficients must be stored as Cnm(n+1,m+1), Snm(n+1,m+1).
%   - The associated Legendre functions here use MATLAB's symbolic
%     legendreP convention. Your coefficient convention must match it.

    arguments
        Cnm (:,:) double
        Snm (:,:) double
        r_ecef (3,1) double
        Re (1,1) double {mustBePositive}
        GM (1,1) double {mustBePositive}
        normalized (1,1) logical = true
    end

    if ~isequal(size(Cnm), size(Snm))
        error('Cnm and Snm must have the same size.');
    end
    if size(Cnm,1) ~= size(Cnm,2)
        error('Cnm and Snm must be square.');
    end

    n_max = size(Cnm,1) - 1;
    n_max = 2;

    % Symbolic variables
    syms x y z real
    rho = sqrt(x^2 + y^2 + z^2);
    s   = z / rho;           % sin(latitude)
    lam = atan2(y, x);

    % Build symbolic potential
    U_sym = sym(0);

    for n = 0:n_max
        radial = (Re / rho)^n;

        inner_sum = sym(0);
        for m = 0:n
            C = Cnm(n+1, m+1);
            S = Snm(n+1, m+1);

            if normalized
                % Fully-normalized associated Legendre function
                Pnm = fullyNormalizedLegendreSym(n, m, s);
            else
                % Unnormalized associated Legendre function
                Pnm = legendreP(n, m, s);
            end

            trig_term = C*cos(m*lam) + S*sin(m*lam);
            inner_sum = inner_sum + Pnm * trig_term;
        end

        U_sym = U_sym + radial * inner_sum;
    end

    U_sym = GM / rho * U_sym;

    % Gradient and Hessian
    a_sym = gradient(U_sym, [x; y; z]);
    T_sym = hessian(U_sym, [x; y; z]);

    % Evaluate at Cartesian point
    x0 = r_ecef(1);
    y0 = r_ecef(2);
    z0 = r_ecef(3);

    U_val = double(vpa(subs(U_sym, [x, y, z], [x0, y0, z0]),50));
    a_val = double(vpa(subs(a_sym, [x, y, z], [x0, y0, z0]),50));
    T_val = double(vpa(subs(T_sym, [x, y, z], [x0, y0, z0]), 50));
end


function Pbar = fullyNormalizedLegendreSym(n, m, s)
% Fully-normalized associated Legendre function built symbolically

    syms t real

    % Ordinary Legendre polynomial P_n(t)
    Pn = legendreP(n, t);

    % Associated Legendre polynomial P_n^m(t)
    % Includes the Condon-Shortley phase (-1)^m
    Pnm = (-1)^m *(1 - t^2)^(m/2) * diff(Pn, t, m);

    % Evaluate at s
    Pnm = subs(Pnm, t, s);

    % Fully-normalized scaling
    delta = double(m == 0);
    Nnm = sqrt((2 - delta) * (2*n + 1) * factorial(sym(n-m)) / factorial(sym(n+m)));

    Pbar = simplify(Nnm * Pnm);
end