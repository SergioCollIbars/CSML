clear;
clc;
close all;
format long g;

addpath('functions/')
%%          CR3BP STM PROPAGATOR
% Description: reconstruct trajectory based on ideal measurements and
% Newton Rapson integatrion method.

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 8*86400;   % [s]
frec = 1/10;        % [Hz]
N = tmax*frec;
t = linspace(tmin, tmax, N);

% define planet parameters
R = 384399e3;  % [m]
m_1 = 5.974E24;  % [Kg]
m_2 = 7.348E22;  % [Kg]
GM = G*(m_1 + m_2); % [m^3 s^-2]
mu =  m_2 / (m_1 + m_2); % mass ratio
planetParams = [mu, R, 1, 1];
poleParams = zeros(6, 1);


% define time dimensionalization
n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
t = t.*n;           % time non-dimensionalization [-]
planetParams(3) = n;

% define initial conditions. L1 orbit
X0 = [1.021968177072928; 0; -0.18206; 0; -0.1031401430288178; 0]; % L1 orbi
% % X0 = [1.0953533743235189E+0; -1.0879975950267760E-28; 0;...
% %      1.3016066486537214E-15; 2.9531900698678965E-1; 0]; % L2 orbit

deltaX0 = [1e-2;0;0;0;0;0];

% % % rotate 2 inertial X0
% % [r0, v0] = rotate2inertial(X0(1:3), X0(4:6), 0, 1);
% % X0 = [r0;v0];

% integatre trajectory and STM
STM_0 = reshape(eye(6,6), [36, 1]);

% define integration options
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% ODE 113
[TIME, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, [], [], "CR3BP"), t, [X0; STM_0], options);
% % [TIME, state_nom] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
% %     poleParams, [], [], "CR3BP_inertial"), t, [X0+deltaX0; STM_0], options);
t = TIME;
N = length(t);


% compute Xdot in CR3BP
Xdot = ones(6, N) * NaN;
for j =1:N
    r = state_true(j, 1:3)';
    v = state_true(j, 4:6)';
    [dU] = computeAcc_CR3BP([r;v]);
    
    Xdot(:, j) = [v', dU];

    
    c = [v', dU]';

    b =  c * [r;v]';
    b =c*c';
    disp(det(b))
end


% % % Reconstruct trajectory from STM
% % for j =1:N
% %     PHI = reshape(state_nom(j, 7:end), [6,6]);
% %     deltaX = PHI * deltaX0;
% %     state_nom(j, 1:6) = state_nom(j, 1:6) - deltaX';
% % end

deltaX0 = (state_true(2, 1:6) - state_true(1, 1:6))';
state_nom = state_true * 0;
state_nom(1, 1:6) = state_true(1, 1:6);
for j =2:N
    PHI = reshape(state_true(j-1, 7:end), [6,6]);
    deltaX = PHI * deltaX0;
    state_nom(j, 1:6) = state_nom(j-1, 1:6) + deltaX';
end

%% PLOT
% plot trajectory
color1 = [204, 0, 204]./256;     % violet
figure()
plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
    'LineWidth', 2, 'Color', 'k')
hold all;
plot3(1-mu, 0, 0, 'LineWidth', 2, 'Marker','o','MarkerFaceColor',color1,'Color',color1)
axis equal
title('3D trajectory in rotating frame')
grid on;
xlabel('X [-]')
ylabel('Y [-]')
zlabel('Z [-]')
legend('S/C', 'Moon')

% % figure()
% % plot3(state_nom(:, 1), state_nom(:, 2), state_nom(:, 3), ...
% %     'LineWidth', 2, 'Color', 'g')
% % axis equal
% % title('Nominal trajectory rotating frame')
% % grid on;
% % xlabel('[-]')
% % ylabel('[-]')
% % zlabel('[-]')


left = [1, 3, 5];
right = [2, 4, 6];
figure()
for j =1:3
    subplot(3, 2, left(j));
    plot(t, state_true(:, j),...
        'LineWidth', 2, 'Color', "#FF00FF")
    hold on;
    plot(t, state_nom(:, j),...
    'LineWidth', 2, 'Color', 'g')

    xlabel('TIME [-]')
    ylabel('r_' + string(j))
    if(j == 1)
        legend('true', 'recons')
    end

    subplot(3, 2, right(j));
    plot(t, state_true(:, j+3), ...
        'LineWidth', 2, 'Color', "#FF00FF")
    hold on;
    plot(t, state_nom(:, j+3), ...
        'LineWidth', 2, 'Color', 'g')
   xlabel('TIME [-]')
   ylabel('v_' + string(j))
end
sgtitle('True trajectory and velocity components. Rotating Frame')

figure()
for j =1:3
    subplot(3, 2, left(j));
    plot(t./86400, abs(state_true(:, j) - state_nom(:, j)),...
        'LineWidth', 2, 'Color', "r")
    xlabel('TIME [-]')
    ylabel('r_' + string(j))
    legend('true - recons')

    subplot(3, 2, right(j));
    plot(t./86400, abs(state_true(:, j+3) - state_nom(:, j+3)), ...
        'LineWidth', 2, 'Color', "r")
   xlabel('TIME [-]')
   ylabel('v_' + string(j))
   legend('true - recons')
end
sgtitle('Trajectory and velocity error components. Rotating Frame')
