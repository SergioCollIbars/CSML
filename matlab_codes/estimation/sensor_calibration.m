clear;
clc;
close all;

%%                     SENSOR CALIBRATION

% design sensor calibration process
addpath('functions/');

% TIME INPUTS
N = 129600;                  % Number of points in the simulation
t_min = 0;                  % Initial time [s]
t_max = 15 * 86400;          % Final time [s]
TIME = linspace(t_min, t_max, N);        % time vector [s]
At = TIME(2) - TIME(1);
Fs = 1/At;

% ORBIT INPUTS
GM = 5.2;           % Gravitational parameter [m^3/s^2]
Re = 246;           % Planet reference radius [m]
a = 2000;           % semi-major axis [m]
e = 0.5;           % eccectricity [-]
i = deg2rad(0);     % inclination [rad]
omega = 0;          % argument of periapsis [rad]
Omega = 0;          % RAAN [rad]
f = 0;              % true anomaly [rad]
rho = a*(1 - e*e);  % orbital parameter [m]

% inital state. IJK frame
r0 = rho / (1 + e * cos(f)) * [cos(f); sin(f); 0];
v0 = sqrt(GM / rho) * [-sin(f); e + cos(f); 0];
[BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
r0 = BN' * r0;
v0 = BN' * v0;
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];


% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM, 0, Re, 1), TIME, X0, options);

r_ACI = state(:, 1:3);
v_ACI = state(:, 4:6);

% construct tensor in RTN + bias + noise
Gamma_RTN = zeros(3*N, 3);
GammaDot_RTN = zeros(3*N, 3);
sigma = 1E-9;     % noise STD [1/s^2]
B = [1, 2, 3;4, 5, 6; 7, 8, 9]*1E-10;

 % flicker noise generation
wn = sigma;
s = 6.32E-12 / wn;      % scaling value
pn = wn / (50 * s);
F = zeros(9, N);
for i = 1:9
    [F(i, :), ~, ~] = noise_profile(wn, pn, N, Fs);
end

% generate measurements. Gamma RTN
Noise = zeros(1, N);
for k = 1:N
    rn = vecnorm(r_ACI(k, :));
    up = 3*k;
    down = up -2;
    E = normrnd(0, sigma, [3, 3]);
    if (k~= 1)
        B = B + ones(3,3)*1E-15 * At;
    end

    Gamma_RTN(down:up, :) = (GM/rn^3)*[-2, 0, 0;0, 1, 0;0, 0, 1] + ...
        E + B + [F(1, k), F(2, k), F(3, k); F(4, k), F(5, k), F(6, k);...
        F(7, k), F(8, k), F(9,k)];

    Noise(k) = E(1, 1) + F(1, 1);

% %     Gamma_RTN(down:up, :) = (GM/rn^3)*[-2, 0, 0;0, 1, 0;0, 0, 1] + ...
% %         E;
% % 
% %     Gamma_RTN(down:up, :) = (GM/rn^3)*[-2, 0, 0;0, 1, 0;0, 0, 1];
end

% generate measurement Gamma Dot RTN
for k = 2:N-1
    up = 3*k;
    down = up -2;

    upP = 3*(k-1);
    downP = upP - 2;
    
    upN = 3*(k+1);
    downN = upN - 2;

    GammaDot_RTN(down:up, :) = (Gamma_RTN(downN:upN, :) - Gamma_RTN(downP:upP, :))./(2*At);
end

% LS estimate biases and drifts. Using Gamma
Ax = 0;
Nx = 0;
R = sigma^2;
Rinv = inv(R);
Y1 = zeros(N, 1);
H1 = zeros(N, 3);
for k = 1:N
    r3 = vecnorm(r_ACI(k, :))^3;
    up =3*k;
    down = up -2;
    g= Gamma_RTN(down:up, :);
    Y1(k) = g(1,1);
    i = k -1;
    y= Y1(k);
    
    H1(k, :) = [-2/r3 ,1, i*At];
    h = H1(k, :);

    % increment Normal equation and information matrix
    Ax = Ax + h' * Rinv * h;
    Nx = Nx + h' * Rinv * y;
end

% solve state
X1 = Ax\Nx;
pos1 = Y1 - H1*X1;

% LS estimate biases and drifts. Using Gamma dot
Ax = 0;
Nx = 0;
R = 100* sigma^2;
Rinv = inv(R);
Y2 = zeros(N, 1);
H2 = zeros(N, 2);
for k = 2:N-2
    rh = r_ACI(k, :)./vecnorm(r_ACI(k, :));
    r4 = vecnorm(r_ACI(k, :))^4;
    v = v_ACI(k, :);
    rv = dot(rh, v);

    up =3*k;
    down = up -2;
    g= GammaDot_RTN(down:up, :);
    Y2(k) = g(1,1);
    y = Y2(k);
    
    H2(k, :) = [6/r4*rv ,1];
    h = H2(k, :);

    % increment Normal equation and information matrix
    Ax = Ax + h' * h;
    Nx = Nx + h' * y;
end

% solve state
X2 = Ax\Nx;
pos2 = Y2 - H2 * X2;

% estimate GM
GM = zeros(1, N);
for k = 1:N-1
    g0= Gamma_RTN(1:3, :);
    Y0 = g0(1,1);
    r03 = vecnorm(r_ACI(1, :))^3;

    up =3*(k);
    down = up -2;
    gk= Gamma_RTN(down:up, :);
    Yk = gk(1,1);
    rk3 = vecnorm(r_ACI(k, :))^3;

    up =3*(k+1);
    down = up -2;
    gk1= Gamma_RTN(down:up, :);
    Yk1 = gk1(1,1);
    rk13 = vecnorm(r_ACI(k+1, :))^3;

    i = k-1;

    A = (Y0 - Yk) / i;
    B = (Y0 - Yk1) / (i+1);

    GM(k) = (A - B)/ (1/(i*r03) - 1/(i*rk3) - 1/(r03*(i+1)) + ...
        1/((i+1)*rk13)) / -2;
end




%%                       PLOT ORBIT

% plot trajectory
figure()
plot3(r_ACI(:, 1), r_ACI(:, 2), r_ACI(:, 3), 'LineWidth', 1.5);
hold on;
plot3(r_ACI(:, 1), r_ACI(:, 2), r_ACI(:, 3), 'o', 'Color', 'r');
axis equal
grid on;
title('Orbit trajectory. IJK')

% plot postfit 1 and 2
figure()
plot(TIME./86400, pos1,'o', TIME./86400, pos2, 'o')
legend('posfit 1', 'postfit 2')
title('Postfit LS')
figure()
plot(TIME./86400, Noise)
title('Noise')

% plot tensor components
G = zeros(9, N);
for j = 1:N
    up = 3 * j;
    down = up -2;
    G(:, j) = reshape(Gamma_RTN(down:up, :), [9, 1]);
end
figure()
for j = 1:9
    subplot(3, 3, j)
    plot(TIME./86400, G(j ,:), 'o', 'Color', 'b')
    xlabel('TIME [days]')
    ylabel('1/s^2')
    grid on;
end

figure()
plot(TIME./86400, GM, '*', 'LineWidth', 1.5)
legend('RMS: ' + string(rms(GM(2:end))))
title('GM estimation value')
%ylim([-6, 6])



