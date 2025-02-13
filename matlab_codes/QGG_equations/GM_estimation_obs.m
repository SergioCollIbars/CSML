clc;
clear;
close all;

% Constants
G = 1; % Gravitational constant (m^3/kg/s^2)
M = 1; % Mass of the point body (kg)

r = 100; % orbit radius
phi = deg2rad(45);  % orbit latitude
w  = deg2rad(1);    % rad/s
t = linspace(0, 1000000, 1000000);

IF_cartesian = 0;
IF_spherical = 0;
IF = zeros(2, length(t));
for j = 1:length(t)
    x = r * sin(phi) * cos(w*t(j));
    y = r * sin(phi) * sin(w*t(j));
    z = r * cos(phi);

    % Variables (x, y, z are coordinates in meters)
    r2 = sqrt(x^2 + y^2 + z^2); % Distance from the point mass to the field point
    
    % Tensor components
    T_xx = G * M / r^3 * (3 * x^2 / r^2 - 1);
    T_yy = G * M / r^3 * (3 * y^2 / r^2 - 1);
    T_zz = G * M / r^3 * (3 * z^2 / r^2 - 1);
    T_xy = G * M / r^3 * (3 * x * y / r^2);
    T_xz = G * M / r^3 * (3 * x * z / r^2);
    T_yz = G * M / r^3 * (3 * y * z / r^2);
    
    % Symmetric tensor
    T = [T_xx, T_xy, T_xz;
         T_xy, T_yy, T_yz;
         T_xz, T_yz, T_zz];
    
    % partials
    At = t(2) - t(1);
    h1 = [T_xx; T_xy; T_xz; T_yy; T_yz; T_zz];
    b2  = eye(6,6);
    d2  = eye(6,6) * At;
    h2 = 1/r^3 * [2;0;0;-1;0;-1];
    
    H1 = [h1, b2, d2];
    H2 = [h2, b2, d2];

    IF_cartesian = IF_cartesian + H1' * H1;
    IF_spherical = IF_spherical + H2' * H2;
    IF(1, j) = det(IF_cartesian);
    IF(2, j) = det(IF_spherical);
end

figure
plot(t, IF, 'LineWidth', 2)
xlabel('TIME [s]')
legend('Cartesian', 'Spherical')
title('Information matrix determinant')

