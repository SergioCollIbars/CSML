clear;
close all;
clc;
format long g;

%%                    SPHERICAL HARMONICS ANALYTIC

% Inputs
rn = 9563E3;                        % orbit radius [m]
phi = 0;                            % latitude [rad]
lambda = pi/2;                      % longitude [rad]

Re = 6378E3;                        % Planet radius [m]   6355E3 (polar radius)
GM = 3.986004418E14;                % gravity param [m^3 / s^2]

n_max = 3;                          % spherical harmonic order

% Hamonics coefficients
C = [1, 0, 0, 0;...
     0, 0, 0, 0;...
     -1.08E-3, 0, 1.57E-6, 0;...
     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7];    % Using Motembruck reference

S = [0, 0, 0, 0;...
     0, 0, 0, 0;...
     0, 0, -9.03E-7, 0;...
     0, 2.68E-7, -2.12E-7, 1.98E-7];         % Using Motembruck reference

% Orbit position
x = rn*cos(phi).*cos(lambda);
y = rn*cos(phi).*sin(lambda);
z = rn*sin(phi);

% Orbit radius vector. Inertial frame
r = [x; y; z];

% spherical Harmonics
option = ["analytical", "cartesian"];
[U_s, dU_s, ddU_s] = sphericalHarmonics(GM, Re, n_max, C, S, ...
        r, option);

disp("potential gradient: ");
disp(dU_s)