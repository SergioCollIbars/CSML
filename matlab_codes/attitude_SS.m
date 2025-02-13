clear;
clc;
close all

%%                     ATTITUDE SIMULATION

% Inlcude folder 
addpath(genpath('./matlab_functions'));

% Input variables

% orbit constants
OMEGA_E = 7.2921159e-5;                       % Earth angular velocity [rad/s]
M_E = 5.972E24;                               % Earth mass [Kg]
G = 6.67E-11;                                 % Gravitational constant [m3 kg-1 s-2]
R_E = 6563e3;                                 % Earth radius [m]

r0 = 1E6 + R_E;                               % orbit radious [m]
p0 = sqrt(G * M_E / (r0^3));                  % orbit angular velocity [rad /s]
   
% Input variables 
t_min = 0;                                    % min time simulation
t_max = 2 * pi / p0;                          % max time simulation
N = 100;                                      % number of points in the simulation

% time vector
t = linspace(t_min, t_max, N);

% compute altitude profile
% generate attitude profile. psi angle
[~, ~, ~, ~, ~, ~, psi, psiDot, psiDdot] = attitude_profile("linear", [0, 0, 3], t, p0);
% generate attitude profile. phi angle
[phi, phiDot, phiDdot, ~, ~, ~, ~, ~, ~] = attitude_profile("sinusoidal", [1, 0, 0], t, 1E-3);
% generate attitude profile. psi2 angle
[~, ~, ~, ~, ~, ~, psi2, psiDot2, psiDdot2] = attitude_profile("sinusoidal", [0, 0, 3], t, 1E-3);

% plot attitude profile
plotProfile(psi, psiDot, psiDdot, phi, phiDot, phiDdot, psi2, psiDot2, psiDdot2, t, '\psi', '\phi', '\psi_2');

% compute angular velocity and angular acc in body axis. 

% Euler transformation
[omega, Omega] = computeAngular_vel_acc(psi, psiDot, psiDdot, phi, phiDot, phiDdot, psi2, psiDot2, psiDdot2, t);

% plot angular vel and acc
plotAng_vel_acc(omega, Omega, t)

% compute the accelartion
a = zeros(3, length(t));
a1 = zeros(3, length(t));
r_inertial = zeros(3, length(t));
r_body = zeros(3, length(t));
for k = 1:length(t)

    A1 = [0.5, 0, 0];  % QGG distance separation x, y, z axis along COM sensor 1 [m]
    A2 = [-0.5, 0, 0]; % QGG distance separation x, y, z axis along COM sensor 1 [m]

    % compute rotational matrix
    [R] = rotMatrix(psi(k), phi(k), psi2(k), [3,1,3]);          % I -> B
    
    % compute plaentary movement
    r_i = r0 * [cos(p0 * t(k)), sin(p0 * t(k)), 0]';   % Inertial frame vector    
    r_b = R * r_i;                                     % Body frame vector
    r_inertial(:, k) = r_i;
    r_body(:, k) = r_b;
    disp(r_i)
    % matricial method. 
    [ad] = compute_acc(G, M_E, omega', Omega', p0 * t(k), vecnorm(r_i), A1, A2, k, R);

    % vectorial method.
    [ad1] = compute_acc_vec(omega, Omega, r_b, A1', A2', k, G, M_E);
    

    a(:, k) = ad;
    a1(:, k) = ad1;

end

% plot acceleration
plotAcc(a, a1, t);

%%                          FUNCTIONS

function [R] = rotMatrix(theta1, theta2, theta3, axis)
    % computes the rotational matrix for a given Euler angles
    
    % matrix definition
    R1 = zeros(3,3);
    R2 = zeros(3,3);
    R3 = zeros(3,3);
    
    % theta 1.
    if(axis(1) == 1)
        R1 = [1, 0, 0;
              0, cos(theta1), sin(theta1);
              0, -sin(theta1), cos(theta1)];
    elseif(axis(1) == 2)
        R1 = [cos(theta1), 0, -sin(theta1);
              0, 1, 0;
              sin(theta1), 0, cos(theta1)];
    elseif(axis(1) == 3)
        R1 = [cos(theta1), sin(theta1), 0;                       
              -sin(theta1), cos(theta1), 0;
              0, 0, 1];
    end

    % theta 2.
    if(axis(2) == 1)
        R2 = [1, 0, 0;
              0, cos(theta2), sin(theta2);
              0, -sin(theta2), cos(theta2)];
    elseif(axis(2) == 2)
        R2 = [cos(theta2), 0, -sin(theta2);
              0, 1, 0;
              sin(theta2), 0, cos(theta2)];
    elseif(axis(2) == 3)
        R2 = [cos(theta2), sin(theta2), 0;                       
              -sin(theta2), cos(theta2), 0;
              0, 0, 1];
    end

    % theta 3.
    if(axis(3) == 1)
        R3 = [1, 0, 0;
              0, cos(theta3), sin(theta3);
              0, -sin(theta3), cos(theta3)];
    elseif(axis(3) == 2)
        R3 = [cos(theta3), 0, -sin(theta3);
              0, 1, 0;
              sin(theta3), 0, cos(theta3)];
    elseif(axis(3) == 3)
        R3 = [cos(theta3), sin(theta3), 0;                       
              -sin(theta3), cos(theta3), 0;
              0, 0, 1];
    end

    % compute R
    R = R3 * R2 * R1;               % WARNING: check expression

end

function [omega, Omega] = computeAngular_vel_acc(theta1, theta1Dot, theta1Ddot, theta2,theta2Dot, theta2Ddot, theta3, theta3Dot, theta3Ddot, t)
    % compute angular velocity and angular acc
    
    % matrix definition
    omega = zeros(3, length(t));
    Omega = zeros(3, length(t));

    for k = 1:length(t)
        % angular velocities
        omega(1, k) = theta1Dot(k) * sin(theta3(k)) * sin(theta2(k)) + cos(theta3(k)) * theta2Dot(k);
        omega(2, k) = cos(theta3(k)) * sin(theta2(k)) * theta1Dot(k) - sin(theta3(k)) * theta2Dot(k);
        omega(3, k) = theta1Dot(k) * cos(theta2(k)) + theta3Dot(k);

        % angular acceleration
        Omega(1, k) = theta1Ddot(k) * (sin(theta3(k)) * sin(theta2(k))) + theta1Dot(k) * (cos(theta3(k)) * sin(theta2(k)) * theta3Dot(k) + sin(theta3(k)) * cos(theta2(k)) * theta2Dot(k))...
                        + theta2Ddot(k) * cos(theta3(k)) - theta2Dot(k) * sin(theta3(k)) * theta3Dot(k);
        Omega(2, k) = theta1Ddot(k) * cos(theta3(k)) * sin(theta2(k)) + theta1Dot(k) * (-sin(theta3(k)) * theta3Dot(k) * sin(theta2(k)) + cos(theta3(k)) * cos(theta2(k)) * theta2Dot(k))...
                        - theta2Ddot(k) * sin(theta3(k)) - theta2Dot(k) * cos(theta3(k)) * theta3Dot(k);
        Omega(3, k) = theta1Ddot(k) * cos(theta2(k)) - theta1Dot(k) * sin(theta2(k)) * theta2Dot(k) + theta3Ddot(k);

    end


end

function plotProfile(theta1, theta1Dot, theta1Ddot, theta2, theta2Dot, theta2Ddot, theta3, theta3Dot, theta3Ddot, t, name1, name2, name3)
    % plot attitude profile
    figure();
    
    subplot(3, 1, 1);
    plot(t, rad2deg(theta1), t, rad2deg(theta2), t, rad2deg(theta3));
    legend(name1, name2, name3);
    title('Euler angles [deg]');
    xlabel('Time [s]');
    ylabel('[\theta_1, \theta_2, \theta_3]');
    grid on;
    
    subplot(3, 1, 2);
    plot(t, rad2deg(theta1Dot), t, rad2deg(theta2Dot), t, rad2deg(theta3Dot));
    legend(name1, name2, name3);
    xlabel('Time [s]');
    ylabel('$[\dot{\theta_1}, \dot{\theta_2}, \dot{\theta_3}]$','interpreter','latex');
    grid on;
    
    subplot(3, 1, 3);
    plot(t, rad2deg(theta1Ddot), t, rad2deg(theta2Ddot), t, rad2deg(theta3Ddot));
    legend(name1, name2, name3);
    xlabel('Time [s]');
    ylabel('$[\ddot{\theta_1}, \ddot{\theta_2}, \ddot{\theta_3}]$','interpreter','latex');
    grid on;
end

function plotAng_vel_acc(omega, Omega, t)
    figure();
    
    subplot(2, 1, 1);
    plot(t, omega(1, :), t, omega(2, :), t, omega(3, :));
    legend('\omega_1', '\omega_2', '\omega_3');
    xlabel('time [s]');
    ylabel('\omega');
    grid on;
    title('Angular vel. Body axis');
    
    subplot(2, 1, 2);
    plot(t, Omega(1, :), t, Omega(2, :), t, Omega(3, :));
    legend('\Omega_1', '\Omega_2', '\Omega_3');
    xlabel('time [s]');
    ylabel('\Omega');
    grid on;
    title('Angular acc. Body axis');
end

function plotAcc(a, a1, t)
figure();

u = mean(vecnorm(a));
scale = 1.1;

plot(t, vecnorm(a))
hold on;
plot(t, vecnorm(a1))
title('Acceletarion module |a|');
xlabel('Time [s]');
ylabel('|a| [m/s^2]');
ylim(scale * [ u - max(vecnorm(a))  , max(vecnorm(a))  + u]);
grid on;

figure();
plot(t, (vecnorm(a) - vecnorm(a1))./ vecnorm(a1) * 100);
title('Relative error');
xlabel('Time [s]');
ylabel('\epsilon %');
grid on;

figure();

u1 = mean(a(1,:));

subplot(3, 1, 1);
title('acceleration components')
plot(t, a(1,:));
xlabel('Time [s]');
ylabel('a_x [m/s^2]');
% ylim(scale * [ u1 - max(a(1,:))  , max(a(1,:))  + u1]);
grid on;

u2 = mean(a(2,:));

subplot(3, 1, 2);
plot(t, a(2,:));
xlabel('Time [s]');
ylabel('a_y [m/s^2]');
% ylim(scale * [ u2 - max(a(2,:))  , max(a(2,:))  + u2]);
grid on;

u3 = mean(a(3,:));

subplot(3, 1, 3);
plot(t, a(3,:));
xlabel('Time [s]');
ylabel('a_z [m/s^2]');
% ylim(scale * [ u3 - max(a(3,:))  , max(a(3,:))  + u3]);
grid on;
end
