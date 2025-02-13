clear;
clc;
close all;

%%                          UKF TEST CODE   
% Description: run the UKF algortihm and test its performance in a simple
% scenario. LEO orbit around Earth using range and range rate measurements.


% Load universe
GM = 3.986004418E14;       % [m^3/s^2]
r0 = [500E3;200E3;500E3];         % [m]
r0n= vecnorm(r0);          % [m]
v0 = [0;sqrt(GM/r0n);0];  % [m/s]
n = sqrt(GM/(r0n^3));      % [rad/s]

T = 2;                  % [rev]
f = 2/1;                % [Hz]
tmin = 0;               % [s]
tmax = T*2*pi/n;        % [s]
Nt = round(tmax * f);   % [-]
t = linspace(tmin, tmax, Nt);

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[time, state_t] = ode113(@(t, x)propagator(t, x, GM), t, [r0;v0], options);

% generate measurements
sigmaRho = 1;               % [m]
sigmaRhoDot = 1E-4;         % [m/s]
R0 = diag([sigmaRho, sigmaRhoDot].^2);
noiseRho = normrnd(0, sigmaRho, [1, Nt]);
noiseRhoDot = normrnd(0, sigmaRhoDot, [1, Nt]);

Y = ones(2, Nt) * NaN;
for j = 1:Nt
    r = state_t(j, 1:3)';
    v = state_t(j, 4:6)';

    Y(:, j) = compute_measuremts(r, v);
end
Y = Y + [noiseRho;noiseRhoDot];

% initial states
err = [10;10;10; 1E-1;1E-1;1E-1];    % [m] and [m/s]
X0 = [r0;v0] + err;
sigmaP = [20;20;20];        % [m]
sigmaV = [.5;.5;.5];           % [m/s]
P0 = diag([sigmaP.^2; sigmaV.^2]);
Ns = 6;
Nm = 2;

% output values
X = zeros(6, Nt); X(:, 1) = X0;
P = zeros(Nt, 6*6); P(1, :) = reshape(P0, [36, 1]);

% UKF routine
for k = 2:Nt
    % create sigma points. State @ k-1
    Xhat_prev = X(:, k-1);
    [xhat_i_prev] = sigmaPoints_state(Ns, P(k-1, :), Xhat_prev);

    % propagate sigma points using N.L funcitons. State @ k
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    time = [t(k-1), t(k)];
    Xhat_i = xhat_i_prev.*0;
    for i = 1:2*Ns   
        [~, state] = ode113(@(t, x)propagator(t, x, GM), time, ...
            xhat_i_prev(:, i), options);
        Xhat_i(:, i) = state(end, 1:6)';
    end

    % apriori state estimate @ time k
    Xhat_min = 1/(2*Ns).* sum(Xhat_i, 2);

    % apriori state covariance @ time k
    A  = zeros(Ns, Ns);
    for i = 1:2*Ns
        A = A + (Xhat_i(:, i) - Xhat_min) * (Xhat_i(:, i) - Xhat_min)';
    end
     P_min = 1/(2*Ns).* A;

     % create sigma points. State @ k
    [Xhat_i] = sigmaPoints_state(Ns, reshape(P_min, [36, 1]), Xhat_min);

    % compute predicted measurements. State @ time k
    Yhat_i = zeros(Nm, 2*Ns);
    for i = 1:2*Ns
        Yhat_i(:, i) = compute_measuremts(Xhat_i(1:3,i), Xhat_i(4:6,i));
    end
    Yhat = 1/(2*Ns).* sum(Yhat_i, 2);

    % compute measurement covariance
     B  = zeros(Nm, Nm);
    for i = 1:2*Ns
        B = B + (Yhat_i(:, i) - Yhat) * (Yhat_i(:, i) - Yhat)';
    end
     Py = 1/(2*Ns).* B + R0;
     
    C  = zeros(Ns, Nm);
    for i = 1:2*Ns
        C = C + (Xhat_i(:, i) - Xhat_min) * (Yhat_i(:, i) - Yhat)';
    end
     Pxy = 1/(2*Ns).* C;

     % kalman update
     K = Pxy/(Py);
     Xhat_plus = Xhat_min + K * (Y(:, k) - Yhat);
     P_plus = P_min - K * Py * K';
    
     % save states
     X(:, k) = Xhat_plus;
     P(k, :) = reshape(P_plus, [36, 1]);
end

sigma = ones(6, Nt);
for j = 1:Nt
    p = reshape(P(j, :), [6,6]);
    sigma(:, j) = sqrt(diag(p));
end

% Plot results
figure()
for j = 1:3
    subplot(3, 1, j)
    plot(t, state_t(:, j)' - X(j, :), 'LineWidth', 2)
    hold on;
    plot(t, sigma(j, :), t, -sigma(j, :), 'LineWidth', 2, ...
         'Color', 'k', 'LineStyle','--');

    xlabel('TIME [sec]')
    ylabel('[m]')
end
sgtitle('Position error along time')

figure()
for j = 1:3
    subplot(3, 1, j)
    plot(t, state_t(:, j+3)' - X(j+3, :), 'LineWidth', 2)
    hold on;
    plot(t, 3*sigma(j+3, :), t, -3*sigma(j+3, :), 'LineWidth', 2, ...
         'Color', 'k', 'LineStyle','--');
    xlabel('TIME [sec]')
    ylabel('[m/s]')
end
sgtitle('Velocity error along time')

% plot orbit
figure()
plot3(state_t(:, 1), state_t(:, 2), state_t(:, 3))
axis equal;
title('Truth trajectory around Earth')

%%                          FUNCTIONS

function [dx] = propagator(t, x, GM)
    r = [x(1);x(2);x(3)];       % [m] 
    a = -GM/(vecnorm(r)^3) * r;  % [m/s]

    dx = [x(4);x(5);x(6);a(1);a(2);a(3)];
end

function [y] = compute_measuremts(r, v)
    % position and velocity vector at time k
    rho = vecnorm(r);
    rhoDot = (r' * v)/rho;
    y = [rho;rhoDot];
end

function [xhat_i] = sigmaPoints_state(Ns, P, Xhat)
    xhat_i = zeros(Ns, 2*Ns);
    xhat_i(:, :) = Xhat.* ones(Ns, 2*Ns);
    % compute matrix square root
    L = chol(Ns.*reshape(P, [6,6]), 'lower');
    X = L';

    % propagate first set of points
    for i = 1:Ns
        xhat_i(:, i)    =  xhat_i(:, i)    + X(i, :)';
        xhat_i(:, i+Ns) =  xhat_i(:, i+Ns) - X(i, :)';
    end
end
