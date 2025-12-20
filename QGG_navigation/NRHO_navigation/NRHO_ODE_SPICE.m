clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")
%%                 QGG ODE INTEGRATOR VS SPICE TEST
% Description: test the ODE 113 integration accuracy vs the SPICE generated
% orbit. Tested orbit: Lunar gateway NRHO. 

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

% Initial configuration
system = "EPHEM";
consider_cov = 0;
tmin = 0;                           % [rad]
tmax = 1*1.4968 + tmin;             % [rad]
frec = 1/60;                        % [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME, DOM] = ...
    load_universe(system, [tmin, tmax], frec);

% load initial conditions & perturbations
X0 = load_initCond(system, planetParams, TIME);
% % deltaX0 = [1E-8, 2E-8, 1E-8, 0,0,0]'; % perturbation 1
% % deltaX0 = [20E-8, 0.1E-8, 1E-8, 0,0,0]'; % perturbation 2
% % deltaX0 = [0, 0, 0, 0,0,0.00001024]'; % perturbation 3
deltaX0 = X0.*[0,0,0,0,0,0]';

% compute time span
n = round((tmax - tmin)/planetParams(3) * frec);
time = linspace(TIME(1), TIME(2), n);

% numerical integration
tol = 1E-13;
for j = 1:length(tol)
    % current tolerance
    tolj = tol(j);

    % integrate trajectory
    options = odeset('RelTol',tolj,'AbsTol',tolj);
    STM0 = reshape(eye(6,6), [36, 1]);
    tStart = tic;
    [t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0, 0), time, ...
        [X0+deltaX0; STM0], options);
    tEnd = toc(tStart);

    % save integrated state
    Xsc_ODE = state(:, 1)';
    Ysc_ODE = state(:, 2)';
    Zsc_ODE = state(:, 3)';

    VXsc_ODE = state(:, 4)';
    VYsc_ODE = state(:, 5)';
    VZsc_ODE = state(:, 6)';
end

TIME = t';   % [-]
jd = 2451545 + TIME / planetParams(3) / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', ...
    'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';

date = humanReadableTime;

% SPICE integration
time = TIME./planetParams(3);    % [s]
[sc_SPICE, ~] = cspice_spkezr('-60000', time, 'J2000', 'NONE', '3');
sc_SPICE(1:3, :) = sc_SPICE(1:3, :)*1E3/planetParams(2); 
sc_SPICE(4:6, :) = sc_SPICE(4:6, :)*1E3/(planetParams(2)*planetParams(3));

[MOON, ~] = cspice_spkezr('MOON', time, 'J2000', 'NONE', '3');
MOON(1:3, :) = MOON(1:3, :)*1E3/planetParams(2); 
MOON(4:6, :) = MOON(4:6, :)*1E3/(planetParams(2)*planetParams(3));

% motion of the primaries
[posE_ODE, posM_ODE] = compute_posPrimaries(TIME, planetParams, system);
cspice_kclear

% plot SPICE trajectory
figure()
plot3(sc_SPICE(1, :) , sc_SPICE(2, :) , sc_SPICE(3, :), 'LineWidth', 2)
hold on;
plot3(Xsc_ODE, Ysc_ODE, Zsc_ODE, 'LineWidth', 2)
axis equal;
grid on;

% plot trajectory components
figure()
subplot(2, 1, 1)
plot(date, vecnorm(sc_SPICE(1:3, :)), 'LineWidth', 2)
title('S/C ID : 6000')

subplot(2, 1, 2)
plot(date, vecnorm([Xsc_ODE;Ysc_ODE;Zsc_ODE]), 'LineWidth', 2)
title('Integrator')

% plot integration error
figure()
subplot(3, 2, 1)
semilogy(date, abs(sc_SPICE(1, :) - Xsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 2, 3)
semilogy(date, abs(sc_SPICE(2, :) - Ysc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 2, 5)
semilogy(date, abs(sc_SPICE(3, :) - Zsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')

subplot(3, 2, 2)
semilogy(date, abs(sc_SPICE(4, :) - VXsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vx [m/s]')
subplot(3, 2, 4)
semilogy(date, abs(sc_SPICE(5, :) - VYsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vy [m/s]')
subplot(3, 2, 6)
semilogy(date, abs(sc_SPICE(6, :) - VZsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vz [m/s]')
xlabel('Epoch [s]')
sgtitle('Component error between SPCIE and filter integrator');


figure()
subplot(3, 1, 1)
plot(date, abs(MOON(1, :) - Xsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 1, 2)
plot(date, abs(MOON(2, :) - Ysc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 1, 3)
plot(date, abs(MOON(3, :) - Zsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')
sgtitle('Position integrated traj. - Moon');


figure()
subplot(3, 1, 1)
plot(date, abs(sc_SPICE(1, :) - MOON(1,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 1, 2)
plot(date, abs(sc_SPICE(2, :) - MOON(2,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 1, 3)
plot(date, abs(sc_SPICE(3, :) - MOON(3,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')
sgtitle('Position SPICE traj. - MOON');


figure()
plot(date, vecnorm(MOON(1:3, :) - [Xsc_ODE;Ysc_ODE;Zsc_ODE])*planetParams(2)/1E3, 'LineWidth', 2)
xlabel('Epoch [s]')
sgtitle('Position integrated traj. - Moon');

figure()
plot(date, vecnorm(sc_SPICE(1:3, :) - MOON(1:3,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('[Km]')
xlabel('date [days]')
sgtitle('Position norm SPICE traj. - MOON');

% save trajectory
traj_pos = [Xsc_ODE;Ysc_ODE;Zsc_ODE]*planetParams(2)/1E3;                        % [km]
traj_vel = [VXsc_ODE;VYsc_ODE;VZsc_ODE]*(planetParams(2)*planetParams(3))/1E3;   % [km/s]
t        = (TIME - TIME(1))./planetParams(3);                                    % [sec]
traj = [traj_pos;traj_vel;t];
name = "perturbation_4_NRHO.mat";
save(name, "traj")

% close SPICE
cspice_kclear