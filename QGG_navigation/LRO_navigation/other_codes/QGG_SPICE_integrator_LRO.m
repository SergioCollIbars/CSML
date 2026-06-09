clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);
addpath('../simplified_functions/');
addpath('../data/');
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data"))
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/functions"))

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_GRAIL.tm')


utc_start = '2012-03-04 00:00:00';
utc_stop  = '2012-03-04 06:00:00';
N         = 4000;           % number of samples
GM_moon   = 4.9028001224453001E12;                       % [m^3/s^2]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

% % tgt       = 'LUNAR RECONNAISSANCE ORBITER';
tgt       = 'GRAIL-B';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

[R_planet] = cspice_bodvrd(observer, 'RADII', 3); % get planet Radius [Km]
R = 1738e3;                            % [m]

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

%% Integrate trajectory with our integrator
[planetParams, Cmat_true, Smat_true] = load_universe();

% initial conditions
X0 = sc_SPICE(1:6, 1);  % [m] & [m/s]

n_maxM = 600;
% % n_maxM = 900;
n_maxE = 2;
% % planetParams(2) = planetParams(2) * (1 - 1.07e-7);

% Integrator
options  = odeset('RelTol',1E-13,'AbsTol',1E-13);
STM0     = reshape(eye(6,6), [36, 1]);
time_vec = [et(1), et(end)]; 
[t, state] = ode113(@(t, x) EOM_LRO_EPHEM_PCT_fixed(t, x, planetParams, ...
    Cmat_true, Smat_true, n_maxM, n_maxE, et(1), et(end)), time_vec, [X0; STM0], options);

[sc_SPICE, ~] = cspice_spkezr(tgt, t', ref_frame, 'NONE', observer);

% convert time
utc  = cspice_et2utc(t', 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

% compute RTN reference frame
BN_mat  = nan(3 * N, 3);
for k = 1:length(t)
    rk = sc_SPICE(1:3, k);
    vk = sc_SPICE(4:6, k);
    [NB] = RTN2ECI(rk, vk);
    
    maxInd = 3 * k; minInd = maxInd - 2;
    BN_mat(minInd:maxInd, :) = NB';
end

error_pos = state(:, 1:3)' - sc_SPICE(1:3, :)*1e3;  % [m]
error_vel = state(:, 4:6)' - sc_SPICE(4:6, :)*1e3;  % [m/s]

figure()
subplot(1, 2, 1)
plot(tUTC, vecnorm(error_pos), 'LineWidth', 2)
xlabel('time'); ylabel('[m]');

subplot(1, 2, 2)
plot(tUTC, vecnorm(error_vel), 'LineWidth', 2)
xlabel('time'); ylabel('[m/s]')

% rotate error in the RNT frame
error_pos_RNT = error_pos.*0; error_vel_RNT = error_vel.*0;
for k = 1:length(t)
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

% compute and plot acceleration error
error_acc_RTN = velocity2acceleration(t, error_vel_RNT);
figure(); tt = ["R", "T", "N"];
for j = 1:3
    subplot(3, 1, j);
    plot(tUTC, error_acc_RTN(j, :), 'LineWidth', 2);
    title(tt(j)); ylabel('[m/s^2]'); grid on;
end
sgtitle('Acceleration error in RTN frame');


% close SPICE
cspice_kclear


%% helpers
function acc = velocity2acceleration(t, vel)
    %VELOCITY2ACCELERATION Compute acceleration from velocity time history.
    %
    % Inputs:
    %   t   : 1 x N or N x 1 time vector [s]
    %   vel : 3 x N velocity history [m/s] or [km/s]
    %
    % Output:
    %   acc : 3 x N acceleration history [same velocity units / s]
    %
    % Uses central differences for interior points and one-sided
    % differences at the boundaries. Allows non-constant dt.
    
    % Ensure time is row vector
    t = t(:).';
    
    % Check dimensions
    if size(vel,2) ~= length(t)
        error('vel must be 3 x N and t must have length N.');
    end
    
    N = length(t);
    acc = nan(size(vel));
    
    if N < 2
        error('At least two time samples are required.');
    end
    
    % Forward difference at first point
    acc(:,1) = (vel(:,2) - vel(:,1)) / (t(2) - t(1));
    
    % Central difference for interior points
    for k = 2:N-1
        acc(:,k) = (vel(:,k+1) - vel(:,k-1)) / (t(k+1) - t(k-1));
    end
    
    % Backward difference at last point
    acc(:,N) = (vel(:,N) - vel(:,N-1)) / (t(N) - t(N-1));
end