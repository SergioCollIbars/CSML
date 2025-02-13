clc;
close all;
clear;
format long g;

%%                      GEOMETRY MATRIX TEST
% Description: This code test the analytical expression for the geometry
% matrix. Given the map coordinates, gravity field order and planet
% parameters the geometry matrix is computed.


% IMPORTS
addpath("functions/");
addpath("//Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/")

% INPUT

% Gravity field
n_max = 2;
Nc = 0;
for j = 1:n_max+1
    Nc  = Nc + j;
end
Nc = Nc - 2;
Ns = Nc - n_max;
Nx = Nc + Ns;
path = "HARMCOEFS_BENNU_CD_1.txt";
[C, S, R, normalized] = readCoeff(path);

% Planet parameters
GM = 5.2;                                                   % [m^3/s^2]

% Map coordinates. ACAF frame
r = 1000;                                                  % [m] 
phi = deg2rad(45);                                         % [rad]
lambda = deg2rad(45);                                      % [rad]
e = 0;                                                     % eccentricity
Ne = 0;                                                    % prime vertical radius
h = r;                                                     % altitude
r_ACAF = [(Ne + h)*cos(phi)*cos(lambda);...
          (Ne + h)*cos(phi)*sin(lambda);...
          (Ne*(1-e^2) + h)*sin(phi)];

% TEST

% transform to ENU coords
ACAF_ENU = [-sin(lambda), -sin(phi)*cos(lambda), cos(phi)*cos(lambda);...
            cos(lambda), -sin(phi)*sin(lambda), cos(phi)*sin(lambda);...
            0, cos(phi), sin(phi)];
ENU_ACAF = ACAF_ENU';

% numerical H value. Visibility matrix
[H_numerical] = potentialGradient_Cnm(n_max, r_ACAF, R, GM, ...
    ENU_ACAF, normalized);

% % H_numerical = [H_numerical(9,:); H_numerical(8,:); H_numerical(7,:);...
% %     H_numerical(5,:); H_numerical(4,:); H_numerical(1,:)];


% analytical value geometric matrix
H_analytical = geometricmatrix(n_max, r, phi, lambda, GM, R, normalized);

Ax_numerical = H_numerical' * H_numerical;
Ax_analytical = H_analytical' * H_analytical;
