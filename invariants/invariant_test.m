clear;
clc;
close all;

%%          INVARIANTS TEST CODE

% paths
addpath('functions/');

% syms
syms e1 e2 e3 e4 e5 e6 e7 e8 e9 t

% Second order gravity tensor. GM values Gamma = f(r, phi, alpha)
Gamma = [2, 1, 4;1, -1, 5;4, 5, -1];
Gamma = [2, 0, 0;0, -1, 0;0, 0, -1] * t;
Gamma = Gamma * 1;

%noise values
eps = [1, 2, 3;4,5,6;7,8,9].*1/10;
eps= [e1, e2, e3;e4, e5, e6;e7, e8, e9];

% measured matrix
M = Gamma + eps;

% rotation matrices
phi = pi/2;
Mx = [1, 0, 0;0, cos(phi),-sin(phi);0,sin(phi), cos(phi)];
My = [cos(phi), 0, sin(phi);0, 1, 0;-sin(phi), 0, cos(phi)];
Mz = [cos(phi), -sin(phi), 0;sin(phi), cos(phi), 0;0, 0, 1];

% rotate matrix
Gammaz = Mz * (Gamma + eps) * Mz';
Gammay = My * (Gamma + eps) * My';
Gammax = Mx * (Gamma + eps) * Mx';

% compute invariants
[I1, I2, I3] = compute_invariant(Gamma+eps);

[~, i2, i3] = compute_invariant(Gammaz);
A1_z = I1;
A2_z = I2 - i2;
A3_z = I3 - i3;
[~, i2, i3] = compute_invariant(Gammay);
A2_y = I2 - i2;
A3_y = I3 - i3;
[~, i2, i3] = compute_invariant(Gammax);
A2_x = I2 - i2;
A3_x = I3 - i3;
B = M(1,2) - M(2,1);
C = M(1, 3) - M(3, 1);
D = M(2,3)- M(3,2);

m = 1000;
a = zeros(2, m);
b = zeros(2, m);
c = zeros(2, m);
d = a;
f = b;
g = c;

s = linspace(0, 2*pi, m);
for j = 1:m
    phi = s(j);
    Mx = [1, 0, 0;0, cos(phi),-sin(phi);0,sin(phi), cos(phi)];
    My = [cos(phi), 0, sin(phi);0, 1, 0;-sin(phi), 0, cos(phi)];
    Mz = [cos(phi), -sin(phi), 0;sin(phi), cos(phi), 0;0, 0, 1];
    
    [~, I2, I3] = compute_invariant(Gamma+eps);

    Gammaz = Mz * (Gamma + eps) * Mz';
    Gammay = My * (Gamma + eps) * My';
    Gammax = Mx * (Gamma + eps) * Mx';
    
    [~, i2, i3] = compute_invariant(Gammaz);
    A2_z = I2 - i2;
    A3_z = I3 - i3;
    [~, i2, i3] = compute_invariant(Gammay);
    A2_y = I2 - i2;
    A3_y = I3 - i3;
    [~, i2, i3] = compute_invariant(Gammax);
    A2_x = I2 - i2;
    A3_x = I3 - i3;

    a(1, j) = A3_x;
    b(1, j) = A3_y;
    c(1, j) = A3_z;

    d(1, j) = A2_x;
    f(1, j) = A2_y;
    g(1, j) = A2_z;
    
    delta = 10;
    [~, I2, I3] = compute_invariant(delta*Gamma+eps);
    
    Gammaz = Mz * (delta*Gamma + eps) * Mz';
    Gammay = My * (delta*Gamma + eps) * My';
    Gammax = Mx * (delta*Gamma + eps) * Mx';
    
    [~, i2, i3] = compute_invariant(Gammaz);
    A2_z = I2 - i2;
    A3_z = I3 - i3;
    [~, i2, i3] = compute_invariant(Gammay);
    A2_y = I2 - i2;
    A3_y = I3 - i3;
    [~, i2, i3] = compute_invariant(Gammax);
    A2_x = I2 - i2;
    A3_x = I3 - i3;

    a(2, j) = A3_x;
    b(2, j) = A3_y;
    c(2, j) = A3_z;

    d(2, j) = A2_x;
    f(2, j) = A2_y;
    g(2, j) = A2_z;
end

figure()
subplot(3, 1, 1)
plot(rad2deg(s), a)
title('\Delta I_3 x component')
grid on;

subplot(3, 1, 2)
plot(rad2deg(s), b)
title('\Delta I_3 y component')
grid on;

subplot(3, 1, 3)
plot(rad2deg(s), c)
title('\Delta I_3 z component')
grid on;

figure()
subplot(3, 1, 1)
plot(rad2deg(s), d)
title('\Delta I_2 x component')
grid on;

subplot(3, 1, 2)
plot(rad2deg(s), f)
title('\Delta I_2 y component')
grid on;

subplot(3, 1, 3)
plot(rad2deg(s), g)
title('\Delta I_2 z component')
grid on;
