clear;
clc;
close all;

A = [1, 2, 3;4,5,6;7,8,9;10,11,12];
B = [4,6,8;2,4,9;1,7,5;9,8,5];

b = [1,2,3]';
a = [4,5,6]';

l = A * a + B * b;
m = length(l);

% Null space B
N  = null(B');

% Projector operator
PN = N * N';

R = diag(([2,3,5,4]).^2);
r1 = inv(N' * R * N);
r2 = inv(PN * R * PN');

% residual with null space
l1 = N' * r1 * l;

% residuals with projector
l2 = PN * r2 * l;

% projected IF matrix
IF1 = (N' * A)' * r1 * (N' * A);
% % IF2 = (PN * A)' * r2 * (PN * A);

Hpr  = B;
hc   = A;
G11 = Hpr' * inv(R) * Hpr;
G22 = hc'  * inv(R) * hc;
G21 = hc'  * inv(R) * Hpr;
G12 = Hpr' * inv(R) * hc;
G2 = G22 - G21 * (G11\G12);

A2 = hc - Hpr * inv(G11) * G12;

PA1 = Hpr * inv(Hpr' * inv(R) * Hpr) * Hpr' * inv(R);

V = eye(m,m) - PA1;
IF3  = (V * A)'  * inv(R) * (V*A);
% %    
