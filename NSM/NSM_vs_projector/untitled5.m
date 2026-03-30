clear;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   NSM vs PROJECTOR: LOSS OF OBSERVABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM = 5.2;
R = 250;

r = 100 + R;
r_vec = r.*[0;0;1];

H = gravityGradient(GM, r_vec(1), r_vec(2), r_vec(3));

% NSM
V = null(H');

% Projector
J = eye(6) - H * ((H'* H) \ H');


%% FUNCTIONS
function H = gravityGradient(GM, x, y, z)
% GRAVITYGRADIENT  Gravity-gradient tensor (Hessian of U = -GM/r).
%   H = gravityGradient(GM, x, y, z)
%   Inputs:
%     GM : gravitational parameter (scalar)
%     x,y,z : position components  (scalars)
%   Output:
%     H : 3x3 gravity-gradient tensor in the same type as inputs

    r2 = x.^2 + y.^2 + z.^2;
    r  = sqrt(r2);
    r5 = r^5;
    f  = GM^2 ./ (r.^3);

    H = f .* [-6/r5*x*y, 6*x*z/r5, 0;...
            -3*(x^2+y^2)/r5, 3*y*z/r5, -3*x*z/r5;...
            -3*y*z/r5, -3*(z^2+x^2)/r5, -3*x*y/r5;...
            6*x*y/r5, 0, -6*y*z/r5;...
            3*x*z/r5, -3*x*y/r5, -3*(z^2+y^2)/r5;...
            0, -6*x*z/r5, 6*y*z/r5];
end