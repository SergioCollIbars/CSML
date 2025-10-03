clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')


utc_start = '2015-05-20 00:00:00';
utc_stop  = '2015-05-30 00:00:00';
N         = 2000;           % number of samples

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'IAU_MOON'; % options: J2000 / IAU_MOON
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
title('LRO Trajectory around the Moon. Moon frame')

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

% close SPICE
cspice_kclear 
