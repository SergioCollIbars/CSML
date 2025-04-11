clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
addpath('../../../QGG_gravEstim/data_files/')
set(0,'defaultAxesFontSize',16);

%%              Null Space Method Geometry analysis
% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 2;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.5E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 5;
f = 1/30;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);
noise0 = zeros(9, Nt);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;reshape(eye(6,6), [36, 1])], options);

GDOP1 = ones(1, Nt) * NaN;
GDOP2 = ones(1, Nt) * NaN;
% compute geometry matrices
 for j = 1:Nt
    % Position and velocity vector (used in NSM)
    rn_ACI = state_t(j, 1:3)';
    vn_ACI = state_t(j, 4:6)';

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % Null space method
    [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI);
    [~, Hc, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ACI', vn_ACI'], ...
            noise0, Cnm, Snm);

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    hp = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)];
    hc =[Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    
    v1 = [hc, hp];
    gdop_trc = diag(inv(v1' * v1)); 
    GDOP1(j) = sum(gdop_trc(1:(Ncs-1)));
    
    v2 = C' * hc;
    gdop_trc = diag(inv(v2'*v2));
    GDOP2(j) = sum(gdop_trc(1:(Ncs-1)));
 end

 figure()
 plot(t, GDOP1, t, GDOP2, 'LineWidth', 2)
 legend('pos+grav', 'NSM')