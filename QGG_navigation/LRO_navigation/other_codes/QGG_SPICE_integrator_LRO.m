clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);
addpath('../simplified_functions/');
addpath('../data/');

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')


utc_start = '2015-03-15 00:00:00';
utc_stop  = '2015-03-15 12:00:00';
N         = 2000;           % number of samples
[GM] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Moon [km^3/s^2]
GM_moon = GM * 1E9;                       % [m^3/s^2]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

[R_planet] = cspice_bodvrd(observer, 'RADII', 3); % get planet Radius [Km]
R = R_planet(1)*1E3;                            % [m]

sc_SPICE(1:3, :) = sc_SPICE(1:3, :).*1E3;    % [m]
sc_SPICE(4:6, :) = sc_SPICE(4:6, :).*1E3;    % [m/s]

figure()
plot3(sc_SPICE(1, :), sc_SPICE(2, :), sc_SPICE(3, :), 'LineWidth', 2)
hold on;

% Create a sphere
[Xs, Ys, Zs] = sphere(100);      % resolution 100x100
surf(R*Xs, R*Ys, R*Zs, ...
    'FaceAlpha', 1, ...          % transparency
    'EdgeColor', 'none', ...
    'FaceColor', [0.8 0.8 0.8]); % light gray

axis equal;
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('LRO Trajectory around the Moon. J2000 frame')

figure()
subplot(1, 2, 1)
plot(et, (vecnorm(sc_SPICE(1:3, :)) - R)./1E3, 'LineWidth', 2)
xlabel('Time');
ylabel('[km]');

subplot(1, 2, 2)
plot(et, vecnorm(sc_SPICE(4:6, :)), 'LineWidth', 2)
xlabel('Time');
ylabel('[m/s]');
sgtitle('postion & velocity norm w.r.t the Moon');

% compute orbital elements
orb_elem = ones(6, N);
for k = 1:N
   [alpha] = orbitalElem(sc_SPICE(1:3, k), sc_SPICE(4:6, k), GM_moon);
   a  = alpha(3);
   n = sqrt(GM / a^3);
   orb_elem(:, k) = [alpha(1), alpha(3)./1E3, wrapTo2Pi(rad2deg(alpha(6:8))),n];
end
figure()
tt = ['e', 'a', 'i', "\omega", "\Omega", 'n'];
time  = (et - et(1))./3600;
for k = 1:6
    subplot(2, 3, k)
    plot(time, orb_elem(k, :), 'LineWidth', 2);
    xlabel('time')
    title(tt(k))
end

% compute RTN reference frame
BN_mat  = nan(3 * N, 3);
for k = 1:N
    rk = sc_SPICE(1:3, k);
    vk = sc_SPICE(4:6, k);
    [NB] = RTN2ECI(rk, vk);
    
    maxInd = 3 * k; minInd = maxInd - 2;
    BN_mat(minInd:maxInd, :) = NB';
end

%% Integrate trajectory with our integrator

[planetParams, Cmat_true, Smat_true] = load_universe();

% initial conditions
X0 = sc_SPICE(1:6, 1);  % [m] & [m/s]

% orbital period (2BP)
% % T = 1/ sqrt(GM / (vecnorm(X0(1:3))^3));

% Integrator
options = odeset('RelTol',1E-13,'AbsTol',1E-13);
STM0 = reshape(eye(6,6), [36, 1]);

[t, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
    Cmat_true, Smat_true), et, [X0; STM0], options);


error_pos = state(:, 1:3)' - sc_SPICE(1:3, :);  % [m]
error_vel = state(:, 4:6)' - sc_SPICE(4:6, :);  % [m/s]

figure()
subplot(1, 2, 1)
plot(et, vecnorm(error_pos), 'LineWidth', 2)
xlabel('time'); ylabel('[m]');

subplot(1, 2, 2)
plot(et, vecnorm(error_vel), 'LineWidth', 2)
xlabel('time'); ylabel('[m/s]')

% rotate error in the RNT frame
error_pos_RNT = error_pos.*0; error_vel_RNT = error_vel.*0;
for k = 1:N
    maxInd = 3 * k; minInd = maxInd - 2;
    BN     =  BN_mat(minInd:maxInd, :);

    error_pos_RNT(:, k) = BN * error_pos(:, k);
    error_vel_RNT(:, k) = BN * error_vel(:, k);
end

figure(); tt = ["R", "T", "N"];
for j = 1:3
    subplot(3, 1, j);
    plot(et, error_pos_RNT(j, :), 'LineWidth', 2);
    title(tt(j)); ylabel('[m]'); grid on;
end
sgtitle('Position error in RTN frame');

figure(); tt = ["R", "T", "N"];
for j = 1:3
    subplot(3, 1, j);
    plot(et, error_vel_RNT(j, :), 'LineWidth', 2);
    title(tt(j)); ylabel('[m/s]'); grid on;
end
sgtitle('Velocity error in RTN frame');

% close SPICE
cspice_kclear
