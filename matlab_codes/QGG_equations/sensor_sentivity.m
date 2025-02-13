clear;
clc;
close all;
%%
% Imports
addpath('./matlab_functions/');
% Study sensor sentivity during J1 estimation. Circular orbit is considered

% Inital conditions. Inertial frame
t_min = 0;
t_max = 6548;
N = 1000;

t = linspace(t_min, t_max, N);

Re = 6378E3;                                    % Earth radious [km]
mu = 3.986E14;

% pos sensors polar coordiantes {r, r_theta, 0}
Ar1 = [0; 0.5; 0];
Ar2 = [0; -0.5; 0];

% Inertial position vector 
x = linspace(Re, 4*Re, N);
y = zeros(1, N);
z = zeros(1, N);

r = [x; y; z];

% Inertial coordinates
X = [1; 0; 0];
Y = [0; 1; 0];

U_i = zeros(1, N);
S_id = zeros(3, N);

for j=1:3
    for k=1:N
        if(j==1)
            Ar1 = [0.5; 0; 0];
            Ar2 = [-0.5; 0; 0];
        elseif(j ==2)
            Ar1 = [0; 0.5; 0];
            Ar2 = [0; -0.5; 0];
        else
            Ar1 = [0; 0; 0.5];
            Ar2 = [0; 0; -0.5];
        end

        % Inertial potential gradient
        U_i(k) = vecnorm(-mu * r(:, k) / ( vecnorm(r(:, k))^3 ));
        
        % polar to inertial
        theta = atan2(r(2, k), r(1, k));
        
        C = [cos(theta), -sin(theta), 0;...
            sin(theta), cos(theta), 0;
            0, 0, 1];
    
        Ar_1 = C' * Ar1;
        Ar_2 = C' * Ar2;
    
        % Inertial position vectors
        r1 = r(:, k) + Ar_1;
        r2 = r(:, k) + Ar_2;
        
        % Sensitivity
        S1 = r1 / (vecnorm(r1)^3);
        S2 = r2 / (vecnorm(r2)^3);
    
        S_id(j, k) = vecnorm(S1 - S2);
    end
end
figure();
plot(vecnorm(r)/1000, U_i);
grid on;
title('SC CoM acceleration in body frame');
ylabel('U_b [m / s^2]');
xlabel('Orbit radious [km]');


figure();
plot(vecnorm(r)/1000, S_id(1, :));
grid on;
title('Distance sensitivity in inertial pointing. Earth');
ylabel('sensitivity [m^{-2}]');
xlabel('Orbit radius [km]');
legend('x-axis Sensor');

figure();
plot(vecnorm(r)/1000, S_id(1, :), vecnorm(r)/1000, S_id(2, :), vecnorm(r)/1000, S_id(3, :));
grid on;
title('Distance sensitivity in inertial pointing. Bennu');
ylabel('sensitivity [m^{-2}]');
xlabel('Orbit radius [km]');
legend('x-axis Sensor', 'y-axis Sensor', 'z-axis Sensor');

figure();
err = S_id(1, :) - S_id(2, :);
plot(vecnorm(r)/1000, err);
grid on;
title('Error between axis sensitivities');
xlabel('Orbit radius [km]');
ylabel('Error [m^{-2}]');

%%                  SENSITIVITY ALONG ORBIT
% orbital parameters
e = 0;                                           % eccentrycity
a = 7563E3;                                      % semi major axis [m]
rho = a * (1 - e^2);                             % orbital param [m]
n = sqrt(mu / a^3);                              % mean motion [rad / s]

i = 0;                                           % inclination [rad]
omega = 0;                                       % arg periapsis [rad]
Omega = 0;                                       % RAAN [rad]
f = 0;                                           % true anomaly [rad]

% compute initial state vectors (orbit frame)
r0 = rho / (1 + e * cos(f)) * [cos(f);...
                              sin(f);...
                              0];
v0 = sqrt(mu / rho) * [-sin(f);...
                      e + cos(f);...
                      0];

% compute rotation matrix
R = rotationMatrix(Omega, i, omega, [3,1,3]);

% rotate initial state vectors. {I, J, K} frame
r0 = R' * r0;
v0 = R' * v0;

% define initial conditions
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];

% define integration options
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

% ODE 113
[~, state] = ode113(@(t, x) EoM(t, x, mu), t, X0, options);
state = state';

% sensitivity
S = zeros(3, length(t));

for k =1:length(t)
    R = rotationMatrix(n * t(k), 0, 0, [3,1,3]);
    
    r = state(1:3, k);

    % polar to inertial
    theta = atan2(r(2), r(1));

    C = [cos(theta), -sin(theta), 0;...
        sin(theta), cos(theta), 0;
        0, 0, 1];

    % X axis
    Ar1 = [0.5;0;0];
    Ar2 = [-0.5;0;0];
 
    Ar1_i = C' * Ar1;
    Ar2_i = C' * Ar2;

    r1 = state(1:3, k) + Ar1_i;
    r2 = state(1:3, k) + Ar2_i;

    S(1, k) = vecnorm(r1/(vecnorm(r1)^3) - r2/(vecnorm(r2)^3));

    % Y axis
    Ar1 = [0;0.5;0];
    Ar2 = [0;-0.5;0];

    Ar1_i = C' * Ar1;
    Ar2_i = C' * Ar2;

    r1 = state(1:3, k) + Ar1_i;
    r2 = state(1:3, k) + Ar2_i;

    S(2, k) = vecnorm(r1/(vecnorm(r1)^3) - r2/(vecnorm(r2)^3));

    % Z axis
    Ar1 = [0;0;0.5];
    Ar2 = [0;0;-0.5];

    Ar1_i = C' * Ar1;
    Ar2_i = C' * Ar2;

    r1 = state(1:3, k) + Ar1_i;
    r2 = state(1:3, k) + Ar2_i;

    S(3, k) = vecnorm(r1/(vecnorm(r1)^3) - r2/(vecnorm(r2)^3));
end

% plots
figure();
plot(t, S(1, :), t, S(2, :), t, S(3, :));
title('Sensitivity along circular orbit');
grid on;
legend('X-axis sensor', 'Y-axis sensor', 'Z-axis sensor');
xlabel('TIME [s]');
ylabel('sensitivity [m^{-2}]');

 %%                     FUNCTIONS
function [dx] = EoM(t, x, mu)
    
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);

    dx = [x(4);...
          x(5);...
          x(6);...
          -mu*x(1)/(r^3);...
          -mu*x(2)/(r^3);...
          -mu*x(3)/(r^3)];
end