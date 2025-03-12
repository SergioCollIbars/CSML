clear;
clc;
close all;
%%                  ATTITUDE SIMULATOR
% Description: Simulates the atittude of the spacecraft in Euler angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SC inertia
I_x = 10;  % kg*m^2
I_y = 15;  % kg*m^2
I_z = 20;  % kg*m^2

I =diag([I_x, I_y, I_z]);

% Initial values nominal
theta0   = zeros(3,1);
thetaDot0 = zeros(3, 1);
state0_nom = [theta0; thetaDot0];

% initial deviations from nominal
delta_theta0 = 1E-7.*[1;0.5;0.7];
delta_thetaDot0 = 1E-10.*[1;1;1];
state0_true = state0_nom + [delta_theta0;delta_thetaDot0];

% time values
time = linspace(0, 100, 100);

% solver
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, state_nom] = ode113(@(t, x) EoM_att(t, x, I, 0), time, state0_nom, options);

% compute new control
M = ones(3, length(time));
for j = 1:length(time)
    psi = state_nom(j , 1); theta = state_nom(j, 2); phi = state_nom(j ,3);
    psiDot = state_nom(j , 4); thetaDot = state_nom(j , 5); phiDot = state_nom(j , 6);
    A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
    omega = A * [psiDot; thetaDot; phiDot];

    M(:, j) = cross(omega, I*omega);
end

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, state_true] = ode113(@(t, x) EoM_att(t, x, I, M), time, state0_true, options);

% plot
figure()
plot(time, (state_true(:, 1:3) - state_nom(:, 1:3))', 'LineWidth', 2)
title('Euler angle error')
ylabel('[rad]');

figure()
plot(time, (state_true(:, 4:6) - state_nom(:, 4:6))', 'LineWidth', 2)
title('Euler angle rate error')
ylabel('[rad/s]');



%% FUNCTION
function [dx] = EoM_att(t, x, I, M)
    % current attitude and attitude rate
    psi = x(1); theta = x(2); phi = x(3);
    psiDot = x(4); thetaDot = x(5); phiDot = x(6);

    % comput angular velocity
     A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
     omega = A * [psiDot; thetaDot; phiDot];

     % compute angular acceleration
     if(M == 0)
        M = cross(omega, I*omega);
     end
     omegaDot = I\(M - cross(omega, I*omega));    % attitude control
     
     % Euler state acceleration
     T22  = -sin(phi)*cos(theta)*phiDot -cos(phi)*sin(theta)*thetaDot;
     T23  = -cos(phi)*cos(theta)*phiDot +sin(phi)*sin(theta)*thetaDot;
     T32  = cos(phi)*sin(theta)*phiDot  +sin(phi)*cos(theta)*thetaDot;
     T33  = -sin(phi)*sin(theta)*phiDot +cos(phi)*cos(theta)*thetaDot;

     T = 1/cos(theta).*[0, sin(phi), cos(phi);...
          0, cos(phi)*cos(theta),-sin(phi)*cos(theta);...
         cos(theta), sin(phi)*sin(theta), cos(phi)*sin(theta)];

     T_dot = cos(theta) * T * sec(theta) * tan(theta) * thetaDot + ...
         1/cos(theta).*[0, cos(phi)*phiDot, -sin(phi)*phiDot;...
                        0, T22, T23;...
                         -sin(theta)*thetaDot, T32, T33];
     thetaDdot  = T_dot * omega + T * omegaDot;

     % time derivative vector
     dx = [psiDot;thetaDot;phiDot;...
         thetaDdot(1);thetaDdot(2);thetaDdot(3);];
end