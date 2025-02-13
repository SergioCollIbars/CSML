clear;
clc;
close all;


%%          INVARIANT TEST 4
% paths
addpath('functions/');

% syms
syms e1 e2 e3 e4 e5 e6 e7 e8 e9 t

% Second order gravity tensor. GM values Gamma = f(r, phi, alpha)
GAMMA = [2, 0, 0;0, -1, 0;0, 0 -1] * t;
Gamma = [6,0, 0;0,-3, 0;0,0,-3]; % t = 3

%noise values
eps = [1, 2, 3;4,5,6;7,8,9];
EPS= [e1, e2, e3;e4, e5, e6;e7, e8, e9];

% rotation matrices
phi = pi/2;
Mx = [1, 0, 0;0, cos(phi),-sin(phi);0,sin(phi), cos(phi)];
My = [cos(phi), 0, sin(phi);0, 1, 0;-sin(phi), 0, cos(phi)];
Mz = [cos(phi), -sin(phi), 0;sin(phi), cos(phi), 0;0, 0, 1];

% omega
omega = Gamma + eps;
OMEGA = GAMMA + EPS;

omega_x = Mx * omega * Mx';
omega_y = My * omega * My';
omega_z = Mz * omega * Mz';

OMEGA_x = Mx * OMEGA * Mx';
OMEGA_y = My * OMEGA * My';
OMEGA_z = Mz * OMEGA * Mz';

[I1_x, I2_x, I3_x] = compute_invariant(omega_x);
[I1_y, I2_y, I3_y] = compute_invariant(omega_y);
[I1_z, I2_z, I3_z] = compute_invariant(omega_z);

[i1_x, i2_x, i3_x] = compute_invariant(OMEGA_x);
[i1_y, i2_y, i3_y] = compute_invariant(OMEGA_y);
[i1_z, i2_z, i3_z] = compute_invariant(OMEGA_z);

% equations
eqn1 = I2_x == i2_x;
eqn2 = I2_y == i2_y;
eqn3 = I2_z == i2_z;
eqn4 = I3_x == i3_x;
eqn5 = I3_y == i3_y;
eqn6 = I3_z == i3_z;
eqn7 = e1 + e5 + e9 == 0;
eqn8 = omega(1, 2) - omega(2, 1) == e2 - e4;
eqn9 = omega(1, 3) - omega(3, 1) == e3 - e7;
eqn10 = omega(2, 3) - omega(3, 2) == e6 - e8;

F = [eqn1; eqn2; eqn3; eqn4; eqn5; eqn6; eqn7; eqn8; eqn9; eqn10];
s = vpasolve(F);
