clear;
clc;
close all;
format long g;

%%               SPACECRAFT GRADIOMETRY CENTER OF MASS SHIFT             %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/29/2024                                                      %
%                                                                         %
%   Description: Code to simulate gradiometer measurements created by the %
%   CoM changes in SC using a spherical shape model.                      %
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
pos2 = [0;0;0];
Nd = 100000;
pos = [linspace(1, 2.50001, Nd);zeros(1, Nd);zeros(1, Nd)];    % [m]
D = 100000;
R = 0.3;              % [m]

mass_SC = ones(1, Nd) * NaN;
U_numerical = ones(1, Nd) * NaN;
U_analytical = ones(1, Nd) * NaN;

Udotmap_numerical = ones(3, Nd) * NaN;
Udotmap_analytical = ones(3, Nd) * NaN;

r_c = ones(3, Nd) * NaN;

for k =1:Nd
    nobj = Nobj + 1;
    dvec = [Dvec; D];
    xvec_0 = [Xvec; pos(1, 1)];
    yvec_0 = [Yvec; pos(2, 1)];
    zvec_0 = [Zvec; pos(3, 1)];

    xvec_f = [Xvec; pos(1, k)];
    yvec_f = [Yvec; pos(2, k)];
    zvec_f = [Zvec; pos(3, k)];

    rvec = [Rvec; R];
    
    [COM_0, mass_0] = compute_COM(nobj, dvec, rvec, xvec_0, yvec_0, zvec_0);
    [COM_f, mass_f] = compute_COM(nobj, dvec, rvec, xvec_f, yvec_f, zvec_f);
    r_c(:, k) = COM_f - COM_0; 

    epsilon = 0.5;
    gradiometerPosA = [0;0;0] + epsilon;
    gradiometerPosB = [0;0;0] - epsilon;
    
    [U0_A, Udot0_A, Uddot0_A] = grav_sphere(nobj, dvec, xvec_0, yvec_0, zvec_0, ...
        rvec, gradiometerPosA);
    [U0_B, Udot0_B, Uddot0_B] = grav_sphere(nobj, dvec, xvec_0, yvec_0, zvec_0, ...
        rvec, gradiometerPosB);

    [Uf_A, Udotf_A, Uddotf_A] = grav_sphere(nobj, dvec, xvec_f, yvec_f, zvec_f, ...
        rvec, gradiometerPosA);
    [Uf_B, Udotf_B, Uddotf_B] = grav_sphere(nobj, dvec, xvec_f, yvec_f, zvec_f, ...
        rvec, gradiometerPosB);
    
    Udotmap_numerical(:, k) = (Udotf_A - Udotf_B) ;
    Udotmap_analytical(:, k) = (Udot0_A - Udot0_B) + (Uddot0_A - Uddot0_B) * r_c(:, k);

    U_numerical(:, k) = (Uf_A - Uf_B) ;
    U_analytical(:, k) =(U0_A - U0_B) -(Udot0_A - Udot0_B)' * r_c(:, k);
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

% CoM change effect
figure()
subplot(1, 3, 1)
plot(vecnorm(r_c),  Udotmap_numerical(1, :), ...
    vecnorm(r_c),  Udotmap_analytical(1, :), 'LineWidth', 2)
xlabel('r_c [m]')
ylabel('\Delta a(R_c + r_c) - \Delta a(R_c) [m/s^2]')

subplot(1, 3, 2)
plot(vecnorm(r_c),  Udotmap_numerical(2, :), ...
    vecnorm(r_c),  Udotmap_analytical(2, :), 'LineWidth', 2)
xlabel('r_c [m]')
ylabel('\Delta a(R_c + r_c) - \Delta a(R_c) [m/s^2]')

subplot(1, 3, 3)
plot(vecnorm(r_c),  Udotmap_numerical(3, :), ...
    vecnorm(r_c),  Udotmap_analytical(3, :), 'LineWidth', 2)
xlabel('r_c [m]')
ylabel('\Delta a(R_c + r_c) - \Delta a(R_c) [m/s^2]')
sgtitle('Acceleration difference between A and B due to CoM shift')

figure()
plot(vecnorm(r_c), U_numerical, ...
    vecnorm(r_c),  U_analytical, 'LineWidth', 2)
xlabel('r_c [m]')
ylabel('\Delta \phi (R_c + r_c) - \Delta \phi (R_c) [m/s^2]')
legend('numerical', 'analytical')


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
