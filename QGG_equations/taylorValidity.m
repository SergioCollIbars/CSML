clc;
close all;
clear;
set(0,'defaultAxesFontSize',16);

%%                  TEST TAYLOR APPROXIMATION VALIDITY

% Description: Testing the approximation of the potential gradient by
% Taylor suring the acc difference computation. All in the inertial fly

% IMPORTS
addpath('functions/');


% INPUTS
GM = 5.2;                               % Gravity param [m^3 s^-2]

r0 = [-1000, 0, 0]';                    % Initial position vector [m];
v0 = [8.83E-18;-0.0721110255092798;0];  % Initial velocity [m];

Nt = 1000;

n = sqrt(GM/  (vecnorm(r0)^3));
T = 2 * pi / n;

t_max = T;
t_min = 0;
t = linspace(t_min, t_max, Nt);

Ax = [1, 0, 0]';

% solve orbit
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM), t, X0, options);
state = state';

% plot orbit
figure();
plot3(state(1, :), state(2, :), state(3, :));
grid on;

% variable save
ad1 = zeros(3, Nt);
ad2 = zeros(3, Nt);

 for k = 1:Nt

    % Inertial coords. CoM
    xi = state(1, k);
    yi = state(2, k);
    zi = state(3, k);
    
    % Inertial coords. Sensor 1
    xi1 = xi + Ax(1)/2;
    yi1 = yi + Ax(2)/2;
    zi1 = zi + Ax(3)/2;

    % Inertial coords. Sensor 2
    xi2 = xi - Ax(1)/2;
    yi2 = yi - Ax(2)/2;
    zi2 = zi - Ax(3)/2;
    
    % potential gradient at each point. Inertial
    U2 = potentialGradient_GM(GM, xi2, yi2, zi2);
    U1  = potentialGradient_GM(GM, xi1, yi1, zi1);

    U1T = potentialGradient2_GM(GM, xi1, yi1, zi1);
    
    % acc difference. Inertial
    ad1(:, k) = (U1 - U2)./Ax(1);

    ad2(:, k) = U1T(1:3, 1);
 end


 err = (ad1 - ad2);
 errR = err./ad1;

 t = t./3600;

 % Plot absolute error
 figure();
  
 plot(t, err, 'LineWidth', 2)
 xlabel(' TIME [h]');
 ylabel('absolute error [1/s^2]');
 grid on;
 legend('\Gamma_{xx}', '\Gamma_{yx}', '\Gamma_{zx}');
 title('Error in first order taylor series');

 % Plot sensor different meas
  figure();

 plot(t, ad1, 'LineWidth', 2)
 xlabel(' TIME [h]');
 ylabel('\Gamma_{ij} [1/s^2]');
 grid on;
 legend('\Gamma_{xx}', '\Gamma_{xy}', '\Gamma_{xz}');
 title('Gradiometer measurements values');