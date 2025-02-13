clear;
clc;
close all;
format long g;

%%                        SPACECRAFT GRADIOMETRY                         %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 12/27/2023                                                      %
%                                                                         %
%   Description: Code to simulate gradiometer measurements created by the %
%   SC gravity using a spherical shape model.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("data/")
addpath("functions/")
set(0,'defaultAxesFontSize',16);

%% CONFIG MODULE
path = "SC_shape_4.txt";

%% SHAPE MODULE
% Read object shape
[Nobj, Xvec, Yvec, Zvec, Rvec, Dvec] = readShape(path);

%% GRADIOMETRY MODULE
% compute gradiometer value at specified position over a grid map
N = 200;
dim = linspace(-6, 6, N);
[Xmap,Ymap] = meshgrid(dim,dim);
Umap = ones(N, N) * NaN;
Udotmap = ones(N, N, 3) * NaN;
Uddotmap = ones(N, N, 9) * NaN;
for j = 1:N % Xmap
    for i = 1:N % Ymap
        gradiometerPos = [Xmap(j, i); Ymap(j, i); 0];
        [U, Udot, Uddot] = grav_sphere(Nobj, Dvec, Xvec, Yvec, Zvec, ...
            Rvec, gradiometerPos);

        Umap(j, i) = U;
        Udotmap(j, i, 1) = Udot(1);
        Udotmap(j, i, 2) = Udot(2);
        Udotmap(j, i, 3) = Udot(3);
        Uddotmap(j, i, :) = reshape(Uddot, [9, 1]);
    end
end

% compute gradiometer value over the three axes
U = ones(1, N, 3) * NaN;
Udot = ones(3, N, 3) * NaN;
Uddot = ones(9, N, 3) * NaN;
for i =1:3
    gradiometerPos = zeros(3, N);
    gradiometerPos(i, :) = dim;
    for j = 1:N
        pos = gradiometerPos(:, j);
        [u, udot, uddot] = grav_sphere(Nobj, Dvec, Xvec, Yvec, Zvec, ...
            Rvec, pos);
        U(1, j, i) = u;
        Udot(:, j, i) = udot;
        Uddot(:, j, i) = reshape(uddot, [9, 1]);
    end
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

% S/C stats
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
disp("S/C total mass = " + string(mass) + " Kg")
disp("S/C total volume = " + string(vol) + " m^3")
disp("S/C CoM = " + string(CoM'))

% Potential map
figure()
surf(Xmap, Ymap, Umap)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('U [J]')
title('S/C gravity potential map')

% Gravity acceleration map
figure();
subplot(1, 3, 1)
surf(Xmap, Ymap, Udotmap(:, :, 1))
xlabel('X [m]')
ylabel('Y [m]')
shading interp
view(0, 90)
title('\nabla U |x')
subplot(1, 3, 2)
surf(Xmap, Ymap, Udotmap(:, :, 2))
shading interp
view(0, 90)
title('\nabla U |y')
subplot(1, 3, 3)
surf(Xmap, Ymap, Udotmap(:, :, 3))
shading interp
view(0, 90)
title('\nabla U |z')
sgtitle('S/C gravity acceleration map')
colorbar

% Gravity SOGT map
figure();
leg = ["xx", "yx","zx", "xy", "yy", "zy", "xz", "yz", "zz"];

for j = 1:9
subplot(3, 3, j)
surf(Xmap, Ymap, Uddotmap(:, :, j))
xlabel('X [m]')
ylabel('Y [m]')
shading interp
view(0, 90)
title('\nabla(\nabla U) ' + leg(j))
colorbar
end
sgtitle('S/C gravity SOGT map')

% potential / acceleration / SOGT over x axis
figure()
A = U(:, :, 1);
B = Udot(:, :, 1);
C = Uddot(:, :, 1);

subplot(2, 1, 1)
plot(dim, A, 'LineWidth', 1.5, 'Color', 'k')
title('U')
xlabel('X [m]')
ylabel('[J]')
subplot(2, 1, 2)
plot(dim, vecnorm(B), 'LineWidth', 1.5, 'Color', 'k')
title('\nabla U')
xlabel('X [m]')
ylabel('[m/s^2]')
sgtitle("Potential and acceleration over X axis")
figure()
for j = 1:9
    subplot(3, 3, j)
    plot(dim, C(j, :), 'LineWidth', 1.5, 'Color', 'k')
    title('\nabla(\nabla U) ' + leg(j))
    xlabel('X [m]')
    ylabel('[1/s^2]')
end
sgtitle("Gradiometer measurements over X axis")

% potential / acceleration / SOGT over y axis
figure()
A = U(:, :, 2);
B = Udot(:, :, 2);
C = Uddot(:, :, 2);

subplot(2, 1, 1)
plot(dim, A, 'LineWidth', 1.5, 'Color', 'k')
title('U')
xlabel('Y [m]')
ylabel('[J]')
subplot(2, 1, 2)
plot(dim, vecnorm(B), 'LineWidth', 1.5, 'Color', 'k')
title('\nabla U')
xlabel('Y [m]')
ylabel('[m/s^2]')
sgtitle("Potential and acceleration over Y axis")
figure()
for j = 1:9
    subplot(3, 3, j)
    plot(dim, C(j, :), 'LineWidth', 1.5, 'Color', 'k')
    title('\nabla(\nabla U) ' + leg(j))
    xlabel('Y [m]')
    ylabel('[1/s^2]')
end
sgtitle("Gradiometer measurements over Y axis")

% potential / acceleration / SOGT over x axis
figure()
A = U(:, :, 3);
B = Udot(:, :, 3);
C = Uddot(:, :, 3);
leg = ["xx", "yx","zx", "xy", "yy", "zy", "xz", "yz", "zz"];
subplot(2, 1, 1)
plot(dim, A, 'LineWidth', 1.5, 'Color', 'k')
title('U')
xlabel('Z [m]')
ylabel('[J]')
subplot(2, 1, 2)
plot(dim, vecnorm(B), 'LineWidth', 1.5, 'Color', 'k')
title('\nabla U')
xlabel('Z [m]')
ylabel('[m/s^2]')
sgtitle("Potential and acceleration over Z axis")
figure()
for j = 1:9
    subplot(3, 3, j)
    plot(dim, C(j, :), 'LineWidth', 1.5, 'Color', 'k')
    title('\nabla(\nabla U) ' + leg(j))
    xlabel('Z [m]')
    ylabel('[1/s^2]')
end
sgtitle("Gradiometer measurements over Z axis")
