clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Null space approach including attitude errors.
% Test assuming inertial nominal attitude control.
% Author: Sergio Coll
% Date: 08/08/25

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

%                Mesurement mask 
%           xx xy xz yx yy yz zx zy zz
mask     =  [1, 1, 1, 0, 1, 1, 0, 0, 1]';

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.36E3;
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
dt = t(2) - t(1);
Nt = length(t);

% attitude error
Amp = 2E-9;
At      = Amp*ones(3, Nt).*[1;0.7;0.5];
dA_dt   = 0.*ones(3, Nt).*[1;0.7;0.5];        % [rad/s]
ddA_ddt = zeros(3, Nt);                         % [rad/s^2]

% attitude nominal value
Amp = 0; frec = 1;
attitude  = Amp.*sin(frec.*t).*ones(3, Nt);                         % nominal attitude [rad]
datt_dt   = Amp.*frec*cos(frec.*t).*ones(3, Nt);
ddatt_ddt = Amp.*-frec^2*sin(frec.*t).*ones(3, Nt); 
[angVel_true, angAcc_true] = compute_angularVals_v2(attitude + At, datt_dt + dA_dt, ddatt_ddt + ddA_ddt);
[angVel_nom, angAcc_nom]   = compute_angularVals_v2(attitude, datt_dt, ddatt_ddt);

% plot S/C attitude
scale = 3600 * 180 / pi;
plot_Attitude(t./T, attitude, At * scale, angVel_nom, ...
    angVel_true - angVel_nom, angAcc_nom, angAcc_true - angAcc_nom);

% noise values from GOCE mission
noise0 = zeros(9, Nt);
sigma1 = 1E-12;
sigma2 = sigma1; sigma3 = sigma1;

means    = zeros(1, 9);
std_devs = [sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3]; 
num_realizations = length(t); % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% contruct ACAF_ACI rotation matrix
ACAF_ACI_mat = zeros(3*Nt, 3) * NaN;
for j = 1:Nt
    Wt = W0 + W * t(j);
    R =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    maxPos = 3 * j; minPos = maxPos - 2;
    ACAF_ACI_mat(minPos:maxPos, :) = R;
end

% output matrices
Yproj      = ones(9, Nt) * NaN;
invariants = ones(2, Nt) * NaN;
for j = 1:Nt
    % RTN rotation matrix
    ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
    % from ACI to Nominal body frame
% %     B_ACI = eye(3,3);
% %     B_ACI = ACI_RTN';

    [ENU_ACAF] = eci2enu_from_r_latlon(ACAF_ACI*rn(:, j));
    B_ACI = ENU_ACAF * ACAF_ACI;
    
    ACAF_B = ACAF_ACI * B_ACI';

    % Null space method correcting for attitude
    [Yc, ~, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);

    T_ACI = reshape(Yc, [3,3])';
    T_BODY = B_ACI * T_ACI * B_ACI';

    Y_BODY = reshape(T_BODY', [9,1]);

    y = Y_BODY(logical(mask));

    [Hrot_grad] = compute_rotPartials_analy(Yc, B_ACI);
% %     [H_num] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, j), ACAF_ACI, ACAF_B);
    Hrot = Hrot_grad(logical(mask), :);
% %     Hrot_num = H_num(logical(mask), :);
    Hrot(abs(Hrot) < 1e-21) = 0;
    
% %     Hrot = zeros(6, 3);
% %     Hrot(3, 2) = GM/(vecnorm(rn(:, j))^3) * 3;
% %     Hrot(2, 1) = -GM/(vecnorm(rn(:, j))^3) * 3;
    
    % null space
    C = null(Hrot');

    Nm = length(C(1, :));
    
    % Compute invariants
    [I1, I2, I3] = compute_invariants(Yc);
    invariants(:, j) = [I2./1E-27;I3./1E-18];
    
    % projected measurement
    Yproj(1:Nm, j) = C' * y;    
end


% plot invariants vs the projected measurements
figure()
tt = ["Y_1", "Y_2", "Y_3"];
time  = t./3600;
for j = 1:3
    subplot(1, 3, j)
    plot(time, Yproj(j, :)./1E-9, 'LineWidth', 2, 'Color', 'g');
    xlabel('Time [hours]')
    ylabel('Eotvos')
    grid on;
    title(tt(j))

end
figure()
tt = ["I_2", "I_3"];

for j = 1:2
    subplot(1, 2, j)
    semilogy(time, invariants(j, :), 'LineWidth', 2, 'Color', 'b');
    xlabel('Time [hours]')
    ylabel('Eotvos^2')
    title(tt(j))
    grid on;
end

%% FUNCTIONS
function [] = plot_Attitude(t, attitude_nom, attitude_true, angVel_nom, ...
    angVel_true, angAcc_nom, angAcc_true)
    figure()
    subplot(3, 2, 1)
    plot(t, attitude_true, 'LineWidth', 2)
    legend('\delta \Psi, yaw', '\delta \theta, pitch', '\delta \phi roll')
    ylabel('Arcsecond')
    title('Error attitude values')
    subplot(3, 2, 3)
    plot(t, rad2deg(angVel_true), 'LineWidth', 2)
    ylabel('[deg/s]')
    legend({'$\delta \omega_3$', '$\delta \omega_2$', '$\delta \omega_1$'}, 'Interpreter', 'latex')
    subplot(3, 2, 5)
    plot(t, rad2deg(angAcc_true), 'LineWidth', 2)
    ylabel('[deg/s^2]')
    legend({'$\delta \dot{\omega_3}$', '$\delta \dot{\omega_2}$', '$\delta \dot{\omega_1}$'}, 'Interpreter', 'latex')
    xlabel('Orbit revolution')
    subplot(3, 2, 2)
    f = wrapTo2Pi(attitude_nom);
    plot(t, rad2deg(f), 'LineWidth', 2)
    legend('\Psi, yaw', '\theta, pitch', '\phi roll')
    ylabel('[deg]')
    title('Nominal attitude values')
    subplot(3, 2, 4)
    plot(t, rad2deg(angVel_nom), 'LineWidth', 2)
    ylabel('[deg/s]')
    legend({'$\omega_3$', '$\omega_2$', '$\omega_1$'}, 'Interpreter', 'latex')
    subplot(3, 2, 6)
    plot(t, rad2deg(angAcc_nom), 'LineWidth', 2)
    ylabel('[deg/s^2]')
    legend({'$\dot{\omega_3}$', '$\dot{\omega_2}$', '$\dot{\omega_1}$'}, 'Interpreter', 'latex')
    xlabel('Orbit revolution')
end

function [I1, I2, I3] = compute_invariants(Y)
    % convert to tensor form
    T = reshape(Y, [3,3])';

    I1  = trace(T);
    I2  = det(T);
    I3  = 0.5 * (trace(T)^2 - trace(T^2));
end

function [R_ENU_ECEF] = eci2enu_from_r_latlon(rE)
% rI   : 3x1 ECI position (meters)
% gmst : Earth rotation angle at epoch (rad), e.g., GMST or ERA
% phi_c: geocentric latitude (rad)
% lambda: longitude (rad, ECEF)
% R_ENU_I: rotation mapping ECI vectors -> local ENU

xE = rE(1); yE = rE(2); zE = rE(3);

% 2) geocentric latitude & longitude
rho = hypot(xE, yE);
phi_c = atan2(zE, rho);
lambda = atan2(yE, xE);

% 3) ENU <- ECEF
sl = sin(lambda); cl = cos(lambda);
sp = sin(phi_c);  cp = cos(phi_c);

R_ENU_ECEF = [ -sl,        cl,       0;
               -sp*cl, -sp*sl,   cp;
                cp*cl,  cp*sl,   sp ];
end
