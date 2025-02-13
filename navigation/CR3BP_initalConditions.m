clear;
clc;
close all;
format long g;

addpath('../../QGG_navigation/functions/integrator/');
addpath('../../QGG_navigation/functions/measurements/');
addpath('../../QGG_navigation/functions/solver/');

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')
%%          CR3BP INITIAL CONDITIONS GENERATOR
% Description: from a given orbit in the CR3BP, generate a new initial
% condition shifted from that initial point. They cna be express in any
% frame: inertial or rotated.

% Initial configuration
system = "CR3BP";    % options: 2BP, CR3BP, F2BP, FCR3BP, EPHEM
phi = deg2rad(0);   % shift angle from X0 vector.

% time parameters
tmin = 0;                              % [rad]
tmax = 1.1*1.3817;                     % [rad]
frec = 1/60;                           % [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME] = ...
    load_universe(system, [tmin, tmax], frec);
mu = planetParams(1);
f_time = frec;                        % fixed integ. frec [Hz]
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);

% load initial conditions. rotating frame (baricenter centered)
r0 = [1.021968177072928; 0; -0.18206];
v0 = [0; -0.1031401430288178; 0]; % L1 orbit

% load initial conditions. inertial frame (baricenter centered)
X0 = load_initCond(system, planetParams);

% basis vectors (recentered)
x = r0(1:3) - [1-mu;0;0];
x_hat = x./vecnorm(x);

% integrate trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
STM0 = reshape(eye(6,6), [36, 1]);

[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0), TIME, [X0; STM0], options);
TIME = t;

state_b = state(:, 1:3)' * 0;
Xb = ones(3, 1) * NaN;
Xi = Xb;
Vi = Xb;
time = NaN;
a = ones(1, length(TIME));
found = 0;
for j = 1:length(TIME)
    % rotation matrix
    NB = [cos(TIME(j)), -sin(TIME(j)), 0;...
        sin(TIME(j)), cos(TIME(j)), 0;...
        0, 0, 1];

    % position vector. Moon 2 S/C. Rotating frame
    r = (NB' * state(j, 1:3)') - [1-mu;0;0];
% %     r = (NB' * state(j, 1:3)');
    r_hat = r./vecnorm(r);
    
    h = cross(x, r);
    h_hat = h./vecnorm(h);

    y_hat = cross(h_hat, x_hat);

    alpha = atan2(dot(r_hat, y_hat), dot(r_hat, x_hat));
    a(j) = alpha;

    if((alpha < (phi + 1.7E-2)) && (alpha > (phi - 1.7E-2)) && (found == 0))
        Xb = r;
        Xi = state(j, 1:3)';
        Vi = state(j, 4:6)';
        time = TIME(j);
        found = 1;
    end

    % rotate to autonomous frame
    state_b(:, j) = r;
end

% print final initial and final points
disp('Initial position in Moon centered rotating frame: ')
disp(x')
disp('Final position in Moon centered rotating frame: ')
disp(Xb')
disp('Final position in barycenter intertial frame: ')
disp(Xi')
disp('Final velocity in barycenter intertial frame: ')
disp(Vi')
disp('Final time in intertial frame: ')
disp(time)

% plot angle vs time
figure()
plot(TIME/planetParams(3)/86400, rad2deg(a), 'LineWidth', 2, 'Color', 'b')
xlabel('Time [days]')
ylabel('\Phi [deg]')
title('Phase angle over time')

% plot in inertial frame. Baricenter
figure()
scale = planetParams(2)/1E3;    % [Km]
plot3(state(:, 1).*scale, state(:, 2).*scale, state(:, 3).*scale, 'LineWidth', 2, 'Color', 'b')
hold all;
plot3(X0(1).*scale, X0(2).*scale, X0(3).*scale, 'Marker','o', 'MarkerFaceColor', 'r' )
plot3(Xi(1).*scale, Xi(2).*scale, Xi(3).*scale, 'Marker','o', 'MarkerFaceColor', 'g')
grid on;
title('Initial and final conditions for NRHO orbit. \Phi = ' + string(phi))
xlabel('X [-]')
ylabel('Y [-]')
zlabel('Z [-]')

% plot in rotating frame. Moon centered
MoonPos = [0;0;0];
arrowInit = [MoonPos, x];
arrowFinal = [MoonPos, Xb];
color1 = '#A020F0';
color2 = 'b';
color3 = 'r';
figure()
plot3(state_b(1, :).*scale, state_b(2, :).*scale, state_b(3, :).*scale, 'LineWidth', 2, 'Color', color1)
hold all;
% % plot3(arrowInit(1, :), arrowInit(2, :), arrowInit(3, :),  'LineWidth', 1.5, 'Color', color2)
% % plot3(arrowFinal(1, :), arrowFinal(2, :), arrowFinal(3, :),  'LineWidth', 1.5, 'Color', color3)
% % plot3(x(1), x(2), x(3),  'Marker', 'o', 'Color', color2, 'MarkerFaceColor', color2)
% % plot3(Xb(1), Xb(2), Xb(3),  'Marker', 'o', 'Color', color3, 'MarkerFaceColor', color3)
plot3(0, 0, 0, 'Marker','o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','k')
% % plot3(-1, 0, 0, 'Marker','o', 'MarkerFaceColor', 'b', 'MarkerSize', 20)
axis equal;
grid on;
legend('', '', '', 'Initial', 'Final', 'Moon', 'Earth')
title('Initial and final conditions for NRHO orbit. \Phi = ' + string(phi))
xlabel('X [Km]')
ylabel('Y [Km]')
zlabel('Z [Km]')
view(85,0)

cspice_kclear