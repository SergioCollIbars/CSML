clear;
clc;
close all;
format long g;

addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);
%%              COMPUTE TENSOR PARTIALS IN THE EIG FRAME
%   Description: compute tensor partials in the eigenvector field using
%   finite differences.


% Harmonic values. A priori
Cnm = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0.1 0.2 0 0 0 0;...
        0.4 0 0 0 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];

% % Cnm = [1 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         -0.4 0 0 0 0 0 0;...
% %         0.4 0 0 0 0 0 0;...
% %         0.02 0 0 0 0 0 0;...
% %         -0.005 0 0 0 0 0 0;...
% %         -0.009 0 0 0 0 0 0];

Cnmap = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0 0 0 0 0 0;...
        0.4 0 0 0 0 0 0;...
        0.02 0 0 0 0 0 0;...
        -0.005 0 0 0 0 0 0;...
        -0.009 0 0 0 0 0 0];

% % Cnmap = [1 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0;...
% %         0 0 0 0 0 0 0];

% % Cnm = Cnmap;
ACnm = Cnm - Cnmap;

Snm = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 1E-6 -3e-05 0 0 0 0;...
        0 0.0001 7e-05 -0.001 0 0 0;
        0 0.007 0.0001 0.0001 0.003 0 0;...
        0 -2e-05 -8e-05 5e-06 -0.0005 0.0005 0;...
        0 -1e-06 -1e-05 -8e-06 -2e-06 -0.0006 -4e-05];

Snmap = zeros(7,7);

% % Snm = Snmap;
ASnm = Snm - Snmap;

% Bennu parameters
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
n_max = 6;
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% Inital nominal conditions
e = 0;                       % eccentrycity
a = 1000;                     % semi major axis [m]
rho = a * (1 - e^2);         % orbital param [m]

i = deg2rad(45);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(45);             % true anomaly [rad]

[r0n, v0n] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);

% Initial true conditions
r0t = r0n + [1;1;1];

% compute nominal GG
Wt = W0 + W * 0;
ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
ACAF_ACI = eye(3,3);
[ACI_RTN] = RTN2ECI(r0n, v0n);
RTN_ACAF = ACI_RTN' * ACAF_ACI';

Cnmref = Cnmap;
[~, ~, ddU_ACAF_nom] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           ACAF_ACI*r0n, Re, GM, 0);

Cnmref(:, 2:end) = Cnmref(:, 2:end).*0;
[~, ~, ddU_ACAF_ref] = potentialGradient_nm(Cnmref, Snmap.*0, n_max, ...
                                           ACAF_ACI*r0n, Re, GM, 0);

% compute meas GG
[~, ~, ddU_ACAF_meas] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           ACAF_ACI*r0t, Re, GM, 0);

EIG_ACAF = RTN_ACAF;
[NB, v] = eig(ddU_ACAF_ref);
BN = NB';
EIG_ACAF = BN;

% compute disturbing GG
[~, ~, ddU_ACAF_ap] = potentialGradient_nm(ACnm, ASnm, n_max, ...
                                           ACAF_ACI*r0n, Re, GM, 0);
[~, ~, ddU_ACAF_ap2] = potentialGradient_nm(ACnm, ASnm, n_max, ...
                                           ACAF_ACI*r0t, Re, GM, 0);

% compute partials
[partials] = computePartials(r0n, Cnm, Snm, GM, Re, n_max, ...
    eye(3,3), ACAF_ACI);

% compute EIG frames
[L2, V2]     = eig(ddU_ACAF_meas);
[L, V]       = eig(ddU_ACAF_ref);
[Lap, Vap]   = eig(ddU_ACAF_ap);

L   = L'; % this gives: BN
L2  = L2';
Lap = Lap';

% compute rotation angles for reference frame
ang1_nom = atan2(L(1, 2), L(1,1));
ang2_nom = -asin(L(1, 3));
ang3_nom = atan2(L(2, 3), L(3,3));

ang1_pos = atan2(L2(1, 2), L2(1,1));
ang2_pos = -asin(L2(1, 3));
ang3_pos = atan2(L2(2, 3), L2(3,3));

Ar = (ang3_pos - ang3_nom); % [rad]
Ap = (ang2_pos - ang2_nom); % [rad]
Al = (ang1_pos - ang1_nom); % [rad]

M3 = [cos(Al), sin(Al), 0;-sin(Al), cos(Al), 0;0, 0, 1];
M2 = [cos(Ap), 0, -sin(Ap); 0, 1, 0;sin(Ap), 0, cos(Ap)];
M1 = [1, 0, 0;0, cos(Ar), sin(Ar);0, -sin(Ar), cos(Ar)];

M = M1 * M2 * M3;

S = M' * V2 * M;    % NOTE: This is the correct way.
D = (V2(3,3) - V2(2,2))* Ar;
j = M' * L2;        % j * ddU_ACAF_meas * j' = S

% project tensor in the EIG frame
ddU_EIG_nom  = EIG_ACAF * ddU_ACAF_nom * EIG_ACAF';
ddU_EIG_ref  = EIG_ACAF * ddU_ACAF_ref * EIG_ACAF';
ddU_EIG_meas = EIG_ACAF * ddU_ACAF_meas * EIG_ACAF';
ddU_EIG_ap   = EIG_ACAF * ddU_ACAF_ap * EIG_ACAF';
ddU_EIG_ap2  = EIG_ACAF * ddU_ACAF_ap2 * EIG_ACAF';

disp(ddU_EIG_meas - ddU_EIG_nom)
disp(partials)

% compute partials at nominal
[Hp] = potentialGradient_Cnm(n_max, ACAF_ACI*r0n, Re, ...
    GM, j, 0);
G3 = V2(3,3);
G2 = V2(2,2);
G1 = V2(1,1);

d = G3*cos(Al)*cos(Ap)*sin(Ap) - G1*cos(Al)*cos(Ap)*sin(Ap);
f = G3*cos(Ap)*sin(Al)*sin(Ap) - G1*cos(Ap)*sin(Al)*sin(Ap);

figure()
L  = L * 50;
L2 = RTN_ACAF;
L2 = L2 * 50;
Lap = Lap * 50;
r0t = ACAF_ACI * r0t;
r0n = ACAF_ACI * r0n;
quiver3(r0n(1),r0n(2),r0n(3),L(1, 3), L(2, 3), L(3,3), 'LineWidth', 2);
hold all;
quiver3(r0n(1),r0n(2),r0n(3),L(1, 2), L(2, 2), L(3,2), 'LineWidth',2);
quiver3(r0n(1),r0n(2),r0n(3),L(1, 1), L(2, 1), L(3,1), 'LineWidth', 2);
legend('nom \lambda 3', 'nom \lambda 2', 'nom \lambda 1')

% % hold all
% % quiver3(r0t(1),r0t(2),r0t(3),L2(1, 3), L2(2, 3), L2(3,3), 'LineWidth', 2);
% % hold all;
% % quiver3(r0t(1),r0t(2),r0t(3),L2(1, 2), L2(2, 2), L2(3,2), 'LineWidth',2);
% % quiver3(r0t(1),r0t(2),r0t(3),L2(1, 1), L2(2, 1), L2(3,1), 'LineWidth', 2);

hold on;
[x,y,z] = sphere;
scale = Re;
x = x * scale;
y = y * scale;
z = z * scale;
surf(x, y, z)
axis equal;

figure()
for k = 1:3
subplot(1, 3, k)
quiver3(r0n(1),r0n(2),r0n(3),L(1, k), L(2, k), L(3, k), 'LineWidth', 2);
hold all;
quiver3(r0n(1),r0n(2),r0n(3),L2(1, k), L2(2, k), L2(3, k), 'LineWidth', 2);
%quiver3(r0n(1),r0n(2),r0n(3),Lap(1, k), Lap(2, k), Lap(3, k), 'LineWidth', 2);
axis equal;
end

function [EIG_ACAF] = computeEIGframe(ddU, dU)
    % compute rotation angle
    th2 = atan2(-2*ddU(1, 2), ddU(1,1)-ddU(2,2));
    th = (th2 -pi)/2;

    R = [cos(th),-sin(th),0;sin(th),cos(th),0;0,0,1];
    
    % compute rotated tensor
    T = R * ddU * R';

    Kmax = -T(1,1)/dU(3);
    Kmin = -T(2,2)/dU(3);
    
    f1 = -T(1, 3)/dU(3);
    f2 = -T(2, 3)/dU(3);

    % compute eigenvalues
    [LL, l] = eig(T);
    
    % compute eigen values
    L = zeros(3,3);
    for k =1:3
        L(:, k) = [-f1*(Kmin - l(k, k));...
                   -f2*(Kmax - l(k,k));...
                   (Kmax-l(k,k))*(Kmin-l(k,k))];
    end
    EIG_ACAF = L * R;
end


function [partials] = computePartials(r0, Cnm, Snm, GM, Re, n_max, ...
    EIG_ACAF, ACAF_ACI)
    % output matrix
    partials = zeros(9, 3);

    % set correctly the coefficient matrix
    Cnm(:, 2:end) = Cnm(:, 2:end).*0; 

    % compute initial spherical coordinates
    r0 = ACAF_ACI * r0; % ACAF coordiantes
    r0n = r0./vecnorm(r0);
    
    phi0 =  atan2(r0n(3), sqrt(r0n(1)^2 + r0n(2)^2));
    lambda0 = atan2(r0n(2), r0n(1));
    rho0 = vecnorm(r0);
    
    for j = 1:3
        Asph = zeros(3, 1);
        % compute displaced spherical coordinates
        Asph(j) = 1E-6;
        rho    = rho0 + Asph(1);
        phi    = phi0 + Asph(2);
        lambda = lambda0 + Asph(3);
        
        rx = rho * cos(phi) * cos(lambda);
        ry = rho * cos(phi) * sin(lambda);
        rz = rho * sin(phi);
        
        rpos = [rx;ry;rz];  % positive increment, ACAF
        rpos = EIG_ACAF * r0 + Asph;
        rpos = EIG_ACAF' * rpos;

        rho    = rho0 - Asph(1);
        phi    = phi0 - Asph(2);
        lambda = lambda0 - Asph(3);
        
        rx = rho * cos(phi) * cos(lambda);
        ry = rho * cos(phi) * sin(lambda);
        rz = rho * sin(phi);
        
        rneg = [rx;ry;rz];  % negative increment, ACAF
        rneg = EIG_ACAF * r0 - Asph;
        rneg = EIG_ACAF' * rneg;

        [~, ~, ddU_ACAF_pos] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rpos, Re, GM, 0);


        [~, ~, ddU_ACAF_neg] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                                   rneg, Re, GM, 0);

        % project into the EIG frame
        ddU_EIG_pos = EIG_ACAF * ddU_ACAF_pos * EIG_ACAF';
        ddU_EIG_neg = EIG_ACAF * ddU_ACAF_neg * EIG_ACAF';

        rpos = EIG_ACAF * rpos;
        rneg = EIG_ACAF * rneg;

        y_pos = reshape(ddU_EIG_pos, [9, 1]);
        y_neg = reshape(ddU_EIG_neg, [9, 1]);
        partials(:, j) = (y_pos - y_neg)./(vecnorm(rpos-rneg));
    end
end