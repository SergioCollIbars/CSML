clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);
addpath('../simplified_functions/');
addpath('../data/');
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data"))
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/functions"))

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')


utc_start = '2012-03-04 00:00:00';
utc_stop  = '2012-03-04 06:00:00';
N         = 4000;           % number of samples
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

% convert time
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

figure()
plot3(sc_SPICE(1, :), sc_SPICE(2, :), sc_SPICE(3, :), 'LineWidth', 2)
hold on; grid on;

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
plot(tUTC, (vecnorm(sc_SPICE(1:3, :)) - R)./1E3, 'LineWidth', 2)
ylabel('[km]');

subplot(1, 2, 2)
plot(tUTC, vecnorm(sc_SPICE(4:6, :)), 'LineWidth', 2)
ylabel('[m/s]');
sgtitle('postion & velocity norm w.r.t the Moon');

% compute orbital elements
orb_elem = ones(6, N);
for k = 1:N
   [alpha] = orbitalElem(sc_SPICE(1:3, k), sc_SPICE(4:6, k), GM_moon);
   a  = alpha(3);
   n = sqrt(GM_moon / (a.^3)); % [rad/s]
   T = (2*pi)./n;
   angles = wrapTo2Pi(alpha(6:8)); % [rad]
   orb_elem(:, k) = [alpha(1), alpha(3)./1E3, rad2deg(angles),T./3600];
end
figure()
tt = ['e', 'a', 'i', "\omega", "\Omega", 'T'];
ylb = ["[-]", "[Km]", "[deg]", "[deg]", "[deg]", "[hr]"];
for k = 1:6
    subplot(2, 3, k)
    plot(tUTC, orb_elem(k, :), 'LineWidth', 2);
    title(tt(k)); ylabel(ylb(k));
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

% Integrator
options  = odeset('RelTol',1E-13,'AbsTol',1E-13);
STM0     = reshape(eye(6,6), [36, 1]);
time_vec = [et(1), et(end)]; 
[t, state] = ode113(@(t, x) EOM_LRO_EPHEM_v0(t, x, planetParams, ...
    Cmat_true, Smat_true), time_vec, [X0; STM0], options);

% interpolate results to SPICE time
x_inter_pos = interp1(t, state(:, 1:3), et', 'spline'); % Nt x 3
x_inter_vel = interp1(t, state(:, 4:6), et', 'spline'); % Nt x 3

error_pos = x_inter_pos' - sc_SPICE(1:3, :);  % [m]
error_vel = x_inter_vel' - sc_SPICE(4:6, :);  % [m/s]

figure()
subplot(1, 2, 1)
plot(tUTC, vecnorm(error_pos), 'LineWidth', 2)
xlabel('time'); ylabel('[m]');

subplot(1, 2, 2)
plot(tUTC, vecnorm(error_vel), 'LineWidth', 2)
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
    plot(tUTC, error_pos_RNT(j, :), 'LineWidth', 2);
    title(tt(j)); ylabel('[m]'); grid on;
end
sgtitle('Position error in RTN frame');

figure(); tt = ["R", "T", "N"];
for j = 1:3
    subplot(3, 1, j);
    plot(tUTC, error_vel_RNT(j, :), 'LineWidth', 2);
    title(tt(j)); ylabel('[m/s]'); grid on;
end
sgtitle('Velocity error in RTN frame');

% close SPICE
cspice_kclear
