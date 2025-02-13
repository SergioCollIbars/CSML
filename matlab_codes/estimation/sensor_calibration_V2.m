clear;
clc;
close all;

%%                     SENSOR CALIBRATION

% design sensor calibration process
addpath('functions/');

% DATA READ
T = readtable('accData.txt');
T2 = readtable('orbitData.txt');
dataAcc = table2array(readtable('accData.txt'));
dataOrbit = table2array(readtable('orbitData.txt'));

TIME = dataAcc(:, 1);
N = length(TIME);
AT = TIME(2) - TIME(1);
G = dataAcc(:, 2:10);
r_ACI = dataOrbit(:, 2:4);
v_ACI = dataOrbit(:, 8:10);

% generate measurements. Gamma RTN
G_RTN = ones(3*N, 3)*NaN;
for k = 1:N
    % rotation matrix
    [ACI_RTN] = RTN2ECI(r_ACI(k, :)', v_ACI(k, :)');
    ACI_RTN = eye(3,3);
    up = 3*k;
    down = up -2;

    % Gamma in ACI frame
    G_ACI = [G(k, 1), G(k, 2), G(k, 3);...
        G(k, 4), G(k, 5), G(k, 6);...
        G(k, 7), G(k, 8), G(k, 9)];
    G_RTN(down:up, :) =  ACI_RTN' * G_ACI * ACI_RTN;
end

% generate measurement Gamma Dot RTN
Gdot_RTN = ones(3*N, 3)*NaN;
for k = 2:N-1
    up = 3*k;
    down = up -2;

    upP = 3*(k-1);
    downP = upP - 2;
    
    upN = 3*(k+1);
    downN = upN - 2;

    Gdot_RTN(down:up, :) = (G_RTN(downN:upN, :) - G_RTN(downP:upP, :))./(2*AT);
end

% LS estimate biases and drifts. Using Gamma
Ax = 0;
Nx = 0;
R = [6.32E-12, 0;0, 2.5E-10];
Rinv = inv(R);
Y1 = zeros(2*N, 1);
H1 = zeros(2*N, 5);
for k = 1:N
    r3 = vecnorm(r_ACI(k, :))^3;
    r5 = vecnorm(r_ACI(k, :))^5;
    x = r_ACI(k, 1);
    y = r_ACI(k, 2);
    up =3*k;
    down = up -2;
    g = G_RTN(down:up, :);
    up = 2*k;
    down = up - 1;
    Y1(down:up) = [g(1,1);g(1, 2)];
    i = k -1;
    yy = Y1(down:up);
    
    H1(down:up, :) = [1/r3-3*x*x/r5 ,1E-8, 1E-20*i*AT, 0, 0;
        1/r3-3*y*y/r5, 0, 0, 1E-8, 1E-20*i*AT];
    h = H1(down:up, :);

    % increment Normal equation and information matrix
    Ax = Ax + h' * Rinv * h;
    Nx = Nx + h' * Rinv * yy;
end

% solve state
X1 = Ax\Nx;
pos1 = zeros(2, N);
for j = 1:N
    up = 2 * j;
    down = up -1;
    pos1(:, j) = Y1(down:up, 1) - H1(down:up, :) * X1;
end
P1 = sqrt(diag(inv(Ax)));

% LS estimate biases and drifts. Using Gamma dot
Ax = 0;
Nx = 0;
sigma = 1E-12;
R = [4.5E-13, 0;0, 4.5E-11].^2;
Rinv = inv(R);
Y2 = zeros(2*N, 1);
H2 = zeros(2*N, 3);
for k = 2:N-2
    rh = r_ACI(k, :)./vecnorm(r_ACI(k, :));
    r4 = vecnorm(r_ACI(k, :))^4;
    r5 = vecnorm(r_ACI(k, :))^5;
    r6 = vecnorm(r_ACI(k, :))^6;
    v = v_ACI(k, :);
    vx = v_ACI(k, 1);
    vy = v_ACI(k, 2);
    x = r_ACI(k, 1);
    y = r_ACI(k, 2);
    rv = dot(rh, v);

    up =3*k;
    down = up -2;
    g= Gdot_RTN(down:up, :);

    up = 2*k;
    down = up - 1;
    Y2(down:up,1) = [g(1,1); g(1, 2)];
    yy = Y2(down:up,1);
    
    H2(down:up, :) = [-3/r4*rv-6*x*vx/r5+15*x*x/r6*rv ,1E-8, 0;...
        -3*(vx*y/r5+vy*x/r5-5*x*y/r6*rv), 0, 1E-8];
    h = H2(down:up, :);

    % increment Normal equation and information matrix
    Ax = Ax + h' * Rinv * h;
    Nx = Nx + h' * Rinv * yy;
end

% solve state
X2 = Ax\Nx;
pos2 = zeros(2, N);
for j = 1:N
    up = 2 * j;
    down = up -1;
    pos2(:, j) = Y2(down:up, 1) - H2(down:up, :) * X2;
end
P2 = sqrt(diag(inv(Ax)));

%% PLOT
figure()
plot(TIME, pos1, 'o', 'Color','b')
title( 'Postfit 1')

figure()
plot(TIME, pos2(1, :), 'o', 'Color','r')
title('Postfit 2, 1')

figure()
plot(TIME, pos2(2, :), 'o', 'Color','r')
title('Postfit 2, 2')


data = [5.19301352256147, 5.19660370925942];
sigma = [0.1, 0.0004];
figure()
errorbar(1:2,data,sigma, 'vertical' ,'LineStyle','none', 'Marker','square', ...
    'LineWidth', 2)
set(gca, 'XTick',1:2, 'XTickLabel',...
            ["Obtained", "Required"])
title('Phase A GM estimation value + \sigma')
grid on;
xlim([-1, 2])
ylabel('[m^3/s^2]')

