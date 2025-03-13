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
theta0     = zeros(3, 1);
thetaDot0  = [0;0;0].*ones(3, 1);
state0_nom = [theta0; thetaDot0];

% initial deviations from nominal
delta_theta0 = 1E-7.*[1;0.5;0.7];
delta_thetaDot0 = 1E-10.*[1;0.5;0.7];
state0_true = state0_nom + [delta_theta0;delta_thetaDot0];

% time values
time = linspace(0, 2*3600, 2*3600);

% solver
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, state_nom] = ode113(@(t, x) EoM_att(t, x, I), time, state0_nom, options);

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, state_true] = ode113(@(t, x) EoM_att(t, x, I), time, state0_true, options);

% plot
figure()
plot(time./3600, (state_true(:, 1:3) - state_nom(:, 1:3))', 'LineWidth', 2);
title('Euler angle - error')
ylabel('[rad]');
grid on;
xlabel('hours')

figure()
plot(time./3600, (state_true(:, 4:6) - state_nom(:, 4:6))', 'LineWidth', 2)
title('Euler angle rate - error')
ylabel('[rad/s]');
grid on;
xlabel('hours')

% compute Euler angle partials
for j = 1:length(time)
    Euler = state_nom(j, 1:3)';
    EulerRate = state_nom(j, 4:6)';
    [dEdtt_dE, dEdtt_dER] = compute_EulerDyn_partials(Euler, EulerRate, I);
end



%% FUNCTION
function [dx] = EoM_att(t, x, I)
    % current attitude and attitude rate
    psi = x(1); theta = x(2); phi = x(3);
    psiDot = x(4); thetaDot = x(5); phiDot = x(6);

    % comput angular velocity
     A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
     omega = A * [psiDot; thetaDot; phiDot];

     % compute angular acceleration
     M = zeros(3, 1);                             % No control
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

     % compute PHI_dot
     % TBD

     % time derivative vector
     dx = [psiDot;thetaDot;phiDot;...
         thetaDdot(1);thetaDdot(2);thetaDdot(3);];
end




function [dEdtt_dE, dEdtt_dER] = compute_EulerDyn_partials(Euler, EulerRate, I)
    eps = 1E-6;
    dEdtt_dE = zeros(3, 3); dEdtt_dER = zeros(3, 3);
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;
        
        % Euler states
        E = Euler + At;
        [E_ddt_pos] = computeEuler_acc(E, EulerRate, I);
        
        E = Euler - At;
        [E_ddt_neg] = computeEuler_acc(E, EulerRate, I);

        % derivative
        dEdtt_dE(:, j) = (E_ddt_pos - E_ddt_neg)./(2.*At(j));
    end

    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;
        
        % Euler states
        ER = EulerRate + At;
        [E_ddt_pos] = computeEuler_acc(Euler, ER, I);
        
        ER = EulerRate - At;
        [E_ddt_neg] = computeEuler_acc(Euler, ER, I);

        % derivative
        dEdtt_dER(:, j) = (E_ddt_pos - E_ddt_neg)./(2.*At(j));
    end
end

function [thetaDdot] = computeEuler_acc(E, EulerRate, I)
    % current attitude and attitude rate
    psi = E(1); theta = E(2); phi = E(3);
    psiDot = EulerRate(1); thetaDot = EulerRate(2); phiDot = EulerRate(3);

    % comput angular velocity
     A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
     omega = A * [psiDot; thetaDot; phiDot];

     % compute angular acceleration
     M = zeros(3, 1);                             % No control
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
end
