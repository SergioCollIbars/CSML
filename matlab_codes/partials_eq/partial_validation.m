clear;
close all;
clc;


%%                       PARTIAL VALIDATION CODE
% ----------------------------------------------------------------------- %
%
% Author: Sergio Coll Ibars
%
% Date: 02/22/2023
%
% Description: Code used to validate partial functions. Get the partial
%   error using Taylor first order.
%
% Inputs: 
%   Math funtion equation and partial equation. 
%   cartesian position and displacement
%
% ----------------------------------------------------------------------- %

% INPUTS
addpath('matlab_functions/');

% INPUTS

% gravity coefficients
GM = 3.98E14;                                               % [m^3 s^-2]
J2 = 1.08E-3;                                               % [??]
Re = 6400E3;                                                % [m]

% nominal coordiantes
x = 4 * Re;                                                 % [m]
y = 4 * Re;                                                 % [m]
z = 4 * Re;                                                 % [m]

% displacement coordinates
deltaX = 1E-5;                                              % [m]
deltaY = 1E-5;                                              % [m]
deltaZ = 1E-5;                                              % [m]


% MAIN
%% GM validation
% compute function at nominal + displacement
dUx = potentialGradient_GM(GM, x+deltaX, y, z);
dUy = potentialGradient_GM(GM, x, y+deltaY, z);
dUz = potentialGradient_GM(GM, x, y, z+deltaZ);

% compute funtion at nominal
dU_n =  potentialGradient_GM(GM, x, y, z);

% compute second partial
ddU = potentialGradient2_GM(GM, x, y, z);

ddU(:, 1) = ddU(:, 1).*deltaX;
ddU(:, 2) = ddU(:, 2).*deltaY;
ddU(:, 3) = ddU(:, 3).*deltaZ;

% compute error matrix
err_GM = [(dUx-dU_n)-ddU(:, 1), (dUy-dU_n)-ddU(:, 2), ...
    (dUz-dU_n)-ddU(:, 3) ];

disp("Error in the GM tensor: ");
disp(err_GM);

%% J2 validation
% compute function at nominal + displacement
dUx = potentialGradient_J2(GM, J2, Re, x+deltaX, y, z);
dUy = potentialGradient_J2(GM, J2, Re, x, y+deltaY, z);
dUz = potentialGradient_J2(GM, J2, Re, x, y, z+deltaZ);

% compute funtion at nominal
dU_n =  potentialGradient_J2(GM, J2, Re, x, y, z);

% compute second partial
ddU = potentialGradient2_J2(GM, J2, Re, x, y, z);

ddU(:, 1) = ddU(:, 1).*deltaX;
ddU(:, 2) = ddU(:, 2).*deltaY;
ddU(:, 3) = ddU(:, 3).*deltaZ;

% compute error matrix
err_J2= [(dUx-dU_n)-ddU(:, 1), (dUy-dU_n)-ddU(:, 2), ...
    (dUz-dU_n)-ddU(:, 3) ];

disp("Error in the J2 tensor: ");
disp(err_J2);


