clear;
clc;
close all;

%%              OBSERVABILITY ANALYSIS FOR NSM NAVIGATION
% Date: 04/18/2025
% Author: Sergio Coll Ibars
% Description: evaluate observability of posiiton and velocity + gravity
% field in asteroid navigation using gradiometer measurements and NSM
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

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
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

% evaluate observability. Information Filter
sigma = 1E-12;
P0 = diag([100, 1E3, 1E3, 1E3, 10, 10, 10].^2);
R0 = diag([sigma, sigma, sigma, sigma, sigma].^2);
PCRB(1, :) = reshape(inv(P0), [1, Nx*Nx]); 
obs = zeros(1, Nt); sigmaP = ones(Nx, Nt) * NaN;
for k = 2:Nt
    PHI_1 = reshape(STM(k-1, :), [Nx, Nx]); % from 0 to k-1
    PHI_2 = reshape(STM(k, :),  [Nx, Nx]);  % from 0 to k
    PHI = PHI_2 * inv(PHI_1);               % from k-1 to k;

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % compute measurements
    [~, Hc, ~] = gradiometer_meas(t(k) ,asterParams, poleParams, [rn(:, k)', vn(:, k)'], ...
                zeros(9, 1), Cnm, Snm, eye(3,3));
    Hc = Hc./GM;
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_ACI);
    
    hc = [Hc(1, :); Hc(4, :); Hc(7, :);Hc(5, :);...
         Hc(8, :)];
    hp = [Hpos(1, :);Hpos(2,:);Hpos(3,:);Hpos(5, :);Hpos(6, :)];
    B = null(hc');
    h = [B' * hc, B' * hp, B' * zeros(5, 3)];
    r = B' * R0 * B;
    
% %     % option 2
% %     h = [hc, hp, zeros(5, 3)];
% %     r = R0;

    Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);
    Ai_min = inv(PHI)' * Aiprev_plus * inv(PHI);
    Ai_plus = Ai_min + h' * (r \ h);
    
    % save observability
    obs(1, k) = rank(Ai_plus);

    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
end

% plot observability
figure()
time = t./3600;
plot(time(2:end), obs(1, 2:end), 'LineWidth', 2, 'Color', 'b')
xlabel('Time [h]')
title('Observability')
grid on; legend('rank', 'condition number')

% compute & propagate state uncertanty
for j = 1:Nt
    p = inv(reshape(PCRB(j, :), [Nx,Nx])); 
    sigmaP(:, j)   = sqrt(diag(p));
end

% plot uncertianty
figure()
subplot(1, 3, 1)
semilogy(time, sigmaP(2:4, :), 'LineWidth', 2)
grid on; ylabel('[m]'); xlabel('Time [h]'); legend('x', 'y', 'z');
title('\sigma position')
subplot(1, 3, 2)
semilogy(time, sigmaP(5:7, :), 'LineWidth', 2)
grid on; ylabel('[m/s]'); xlabel('Time [h]'); legend('x', 'y', 'z');
title('\sigma velocity')
subplot(1, 3, 3)
semilogy(time, sigmaP(1, :), 'LineWidth', 2)
grid on; ylabel('[m^3/s^2]'); xlabel('Time [h]');
title('\sigma GM')