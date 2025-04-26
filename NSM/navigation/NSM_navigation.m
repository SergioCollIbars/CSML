clear;
clc;
close all;

%%              NSM NAVIGATION
% Date: 04/18/2025
% Author: Sergio Coll Ibars
% Description: try navigation using the EKF algorith + the NSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 0;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);
Nx = 6 + Ncs;

% Initial conditions
r      = 1E3;
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
rev = 1;
f = 1/10;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6 + Ncs, 6 + Ncs), [(6+Ncs)^2, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [r0;v0;PHI0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';
STM = state_t(:, 7:end);

% generate noise
sigma = 1E-12;
noise = normrnd(0, sigma, [9, Nt]);
R0 = diag([sigma, sigma, sigma, sigma, sigma].^2);

% generate measurements
[Ytrue, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, [rn', vn'], ...
                noise, Cnm, Snm, eye(3,3));

sigma_n  = 2;
Xp       = abs(normrnd(0, sigma_n, [1,1]));
Cp = Cnm;         Sp = Snm;                                                                % apriori uncertanty gra. field
sigmaPos = [1E1;1E1;1E1;1E-1;1E-1;1E-1];                                                % apriori uncertainty s/c state
P0 = sigma_n^2;

% evaluate observability. Information Filter
P = blkdiag(P0, diag(sigmaPos.^2));
s = state_t(1, 1:6)';
for k = 2:Nt
    asterParams = [Xp, Re, n_max, normalized];

     % time span
    t_span = [t(k - 1), t(k)];

    % integrate trajectory
     [~, state_n] = ode113(@(t, x) EoM(t, x, Cp, Sp, n_max, Xp, Re, normalized, ...
    W0, W, RA, DEC, 1), t_span, [s;reshape(eye(Nx,Nx), [Nx*Nx, 1])], options);
    STM = state_n(:, 7:end);

    % current states 
    rn = state_n(end, 1:3);
    vn = state_n(end, 4:6);

    PHIt = reshape(STM(end, :), [Nx,Nx]);

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % EKF method
    [~, Hc, ~] = gradiometer_meas(t(k) ,asterParams, poleParams, [rn(:, k)', vn(:, k)'], ...
                zeros(9, 1), Cnm, Snm, eye(3,3));
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, Xp, rn(:, k), ACAF_ACI, ACAF_ACI);
    
    hc = [Hc(1, :); Hc(4, :); Hc(7, :);Hc(5, :);...
         Hc(8, :)];
    hp = [Hpos(1, :);Hpos(2,:);Hpos(3,:);Hpos(5, :);Hpos(6, :)];
    B = null(hc');

    dY = B'*(Ytrue(:, k) - Yc);
    % FINISH ... 
    [Xhat, P] = EKF_method(dY, HC_ACI, Hpos, R_N, PHIt, noise(:, j), P);
    s = [rn,vn]' + Xhat(46:end);
    Xp_E(2:end) = Xp_E(2:end) + Xhat;
    [Cp_E, Sp_E] = list2mat(n_max, Nc, Ns, Xp_E(1:46));
    PHI_1 = reshape(STM(k-1, :), [Nx, Nx]); % from 0 to k-1
    PHI_2 = reshape(STM(k, :),  [Nx, Nx]);  % from 0 to k
    PHI = PHI_2 * inv(PHI_1);               % from k-1 to k;

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % compute measurements

    h = [B' * hc, B' * hp, B' * zeros(5, 3)];
    r = B' * R0 * B;
    
    
    Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);
    Ai_min = inv(PHI)' * Aiprev_plus * inv(PHI);
    Ai_plus = Ai_min + h' * (r \ h);
    
    % save observability
    obs(1, k) = rank(Ai_plus);

    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
end