clear;
clc;
close all;
format long g;

%%               SPACECRAFT GRADIOMETRY MASS LOSS EFFECT                 %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/29/2024                                                      %
%                                                                         %
%   Description: Code to simulate gradiometer measurements created by the %
%   mass changes in SC using a spherical shape model.                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("data/")
addpath("functions/")
set(0,'defaultAxesFontSize',16);

%% CONFIG MODULE
path = "SC_shape_3.txt";

%% SHAPE MODULE
% Read object shape
[Nobj, Xvec, Yvec, Zvec, Rvec, Dvec] = readShape(path);

%% GRADIOMETRY MODULE
% compute gradiometer value at specified position over a grid map

% Extra mass definition
pos = [0.5;0.5;0];    % [m]
pos2 = [0.2;0.2;0.2];
Nd = 100;
D = linspace(50, 0, Nd);
R = 0.3;              % [m]

mass_SC = ones(1, Nd) * NaN;
Udotmap = ones(3, Nd) * NaN;

for k =1:Nd     
    nobj = Nobj + 1;
    dvec = [Dvec; D(k)];
    xvec = [Xvec; pos(1)];
    yvec = [Yvec; pos(2)];
    zvec = [Zvec; pos(3)];
    rvec = [Rvec; R];

    [COM, mass_SC(k)] = compute_COM(nobj, dvec, rvec, xvec, yvec, zvec);
    epsilon = 0.5;
    gradiometerPos1 = pos2 + epsilon;
    gradiometerPos2 = pos2 - epsilon;
    
    [~, Udot1, ~] = grav_sphere(nobj, dvec, xvec, yvec, zvec, ...
        rvec, gradiometerPos1);

    [~, Udot2, ~] = grav_sphere(nobj, dvec, xvec, yvec, zvec, ...
        rvec, gradiometerPos2);
    
    Udotmap(1, k) = Udot1(1) - Udot2(1);
    Udotmap(2, k) = Udot1(2) - Udot2(2);
    Udotmap(3, k) = Udot1(3) - Udot2(3);
end


%% POSTPROCESS
% Plot spacecraft shape
figure();
for n = 1:Nobj
    [x,y,z] = sphere;
    
    X = x * Rvec(n) + Xvec(n);
    Y = y * Rvec(n) + Yvec(n);
    Z = z * Rvec(n) + Zvec(n);
    surf(X, Y, Z);
    axis equal
    hold on;
end
title('Spacecraft (S/C) shape. Sphere model')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

% density change effect
figure()
subplot(1, 3, 1)
plot(D,  Udotmap(1, :) -  Udotmap(1, end), 'LineWidth', 2)
xlabel('Density \rho [kg/m^3]')
ylabel('\Delta a_A - \Delta a_B [m/s^2]')

subplot(1, 3, 2)
plot(D,  Udotmap(2, :) -  Udotmap(2, end), 'LineWidth', 2)
xlabel('Density \rho [kg/m^3]')
ylabel('\Delta a_A - \Delta a_B [m/s^2]')

subplot(1, 3, 3)
plot(D,  Udotmap(3, :) -  Udotmap(3, end), 'LineWidth', 2)
xlabel('Density \rho [kg/m^3]')
ylabel('\Delta a_A - \Delta a_B [m/s^2]')
sgtitle('Acceleration difference between A and B due to density change')

% S/C mass change
figure()
plot(D, mass_SC, 'LineWidth', 2)
ylabel('S/C  [kg]')
xlabel('S/C deposit change')
title('Total mass of the S/C')

%% FUNCTIONS
function [CoM, mass] = compute_COM(Nobj, Dvec, Rvec, Xvec, Yvec, Zvec)
    mass = 0;
    vol = 0;
    CoM_x = 0;
    CoM_y = 0;
    CoM_z = 0;
    for n = 1:Nobj
        mass = mass + Dvec(n) * 4/3 * pi * Rvec(n)^3;
        vol = vol + 4/3 * pi * Rvec(n)^3;
        CoM_x = CoM_x + Dvec(n) * 4/3 * pi * Rvec(n)^3 * Xvec(n);
        CoM_y = CoM_y + Dvec(n) * 4/3 * pi * Rvec(n)^3 * Yvec(n);
        CoM_z = CoM_z + Dvec(n) * 4/3 * pi * Rvec(n)^3 * Zvec(n);
    end
    CoM = [CoM_x; CoM_y; CoM_z]./mass;
end
