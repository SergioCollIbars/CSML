clear;
clc;
close all;
addpath('functions/');

%%          PROVE INVARIANTS UNDER FRAME ROTATION FOR GRAV. TENSOR

% gravity tensor
GM = 3.986004418E14;
R = 6.3781E6;
r = 2000e3 + R;
t = GM/(r^3);

Gamma  = [2, 0, 0;0, -1, 0;0, 0, -1] * t;
Gamma2 = [4, 5, 6;5, -3, 3;6, 3, -1] * t;
dGamma = [-6, 0, 0;0, 3, 0;0, 0, 3] * t/r;

delta  = 1020;
r = r + delta;
t = GM/(r^3);
deltaGamma = [2, 0, 0;0, -1, 0;0, 0, -1] * t;

% rotation matrix
phi = pi/4;
Mx = [1, 0, 0;0, cos(phi),-sin(phi);0,sin(phi), cos(phi)];
BN = Mx;

% invariants
[I1,I2, I3] = compute_invariant(BN*Gamma*BN');
[I1n,I2n, I3n] = compute_invariant(Gamma);

% partials in RTN frame
dI2_dr = dGamma(1,1)*(Gamma(2,2) + Gamma(3,3)) + ...
    Gamma(1,1) * (dGamma(2,2) + dGamma(3,3))   + ...
    dGamma(2,2)*Gamma(3,3) + Gamma(2,2)*dGamma(3,3);

A = I2 - I2n;
B = dI2_dr;
X = inv(B) * A;


% sum of tensors
dU = Gamma + Gamma2;
[I1,I2, I3] = compute_invariant(dU);

[I1n,I2n, I3n] = compute_invariant(Gamma2);
[I1n2,I2n2, I3n2] = compute_invariant(Gamma2 + Gamma);

D = I2n2 + I2n;