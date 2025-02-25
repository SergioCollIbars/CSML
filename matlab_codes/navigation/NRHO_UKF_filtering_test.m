clear;
clc;
close all;

%%                          UKF TEST IN NRHO CODE   
% Description: run the UKF algortihm and test its performance the NRHO orbit 
% using gradiometer measurements.


% Load universe
D = 384399e3;               % [m]
mu = 0.0121505843958292;    % [-] 
n =  2.66532477988268e-06;  % [1/s]

r0 = [1.021968177072928; 0; -0.18206];
v0 = [0; -0.1031401430288178; 0]; % L1 orbit
P = 1.4968;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = load('EM_NHalo_L2_Family.mat');
j = 50;
x0   = data.trajFam{j,1}.iState;
r0 = x0(1:3); v0 = x0(4:6);
mu = data.trajFam{j,1}.mu;
P    = data.trajFam{j, 1}.period;

f = 1/60;                   % [Hz]
tmin = 0;                   % [-]
tmax = 3*P;                 % [-]
Nt = round(tmax * f / n);   % [-]
time = linspace(tmin, tmax, Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r0, v0] = rotate2inertial(r0, v0, 0, 1);  % [-] and [-]

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x)propagator(t, x, mu), time, [r0;v0], options);

% generate measurements
sigma = 1E-12;                       % [1/s^2]
sigmaVec = ones(6, 1).*sigma/(n^2);  % [-]
R0 = diag(sigmaVec.^2);
noise = normrnd(0, sigma, [6, Nt]);

Y = ones(6, Nt) * NaN;
for j = 1:Nt
    Y(:, j) = compute_measuremts(time(j), state_t(j, 1:6)', mu, ...
        noise(:, j));
end

% initial states
errP = 0.*[50;50;50]./D;                   % [-]
errV = 0.*[1E-1;1E-1;1E-1]./(D*n);         % [-]
X0 = [r0;v0] + [errP;errV];
sigmaP = [100;100;100]./D;                 % [-]
sigmaV = [.5;.5;.5]./(D*n);             % [-]
P0 = diag([sigmaP.^2; sigmaV.^2]);
Ns = 6;
Nm = 6;

% output values
X = zeros(6, Nt); X(:, 1) = X0;
P = zeros(Nt, 6*6); P(1, :) = reshape(P0, [36, 1]);

% UKF routine
for k = 2:Nt
    disp('Loading ...' + string(k/Nt*100) + '%')
    % create sigma points. State @ k-1
    Xhat_prev = X(:, k-1);
    [xhat_i_prev] = sigmaPoints_state(Ns, P(k-1, :), Xhat_prev);

    % propagate sigma points using N.L funcitons. State @ k
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    t = [time(k-1), time(k)];
    Xhat_i = xhat_i_prev.*0;
    for i = 1:2*Ns   
        [~, state] = ode113(@(t, x)propagator(t, x, mu), t, ...
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
        Yhat_i(:, i) = compute_measuremts(t(end), Xhat_i(1:6,i), mu, zeros(6, 1));
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
     % % Xhat_plus = Xhat_min + K * (Y(:, k) - Yhat);
     Xhat_plus = Xhat_min;
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
    semilogy(time, abs(state_t(:, j)' - X(j, :)), 'LineWidth', 2)
    hold on;
    semilogy(time, sigma(j, :),  'LineWidth', 2, ...
         'Color', 'k', 'LineStyle','--');

    xlabel('TIME [sec]')
    ylabel('[m]')
end
sgtitle('Position error along time')

figure()
for j = 1:3
    subplot(3, 1, j)
    semilogy(time, abs(state_t(:, j+3)' - X(j+3, :)), 'LineWidth', 2)
    hold on;
    semilogy(time, 3*sigma(j+3, :), 'LineWidth', 2, ...
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

function [dx] = propagator(t, x, mu)
        M = t;
        r1 = [x(1)+mu*cos(M);x(2)+mu*sin(M);x(3)];              % SC-Earth
        r2 = [x(1)+(mu - 1)*cos(M); x(2)+(mu-1)*sin(M); x(3)];  % SC-Moon
        
        Cmat = zeros(7, 7); Smat = Cmat;
        Cmat(1, 1) = 1; Smat(1, 1) = 0;
        n_max = 0; normalized = 1; Re = 1;

        [~, dU1, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r1, Re, 1-mu, ...
                                                    normalized);

        [~, dU2, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r2, Re, mu, ...
                                                    normalized);
        dU = (dU1 + dU2);

        dx = [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3)];
end

function [y] = compute_measuremts(t, x, mu, noise)
        M = t;
        r1 = [x(1)+mu*cos(M);x(2)+mu*sin(M);x(3)];              % SC-Earth
        r2 = [x(1)+(mu - 1)*cos(M); x(2)+(mu-1)*sin(M); x(3)];  % SC-Moon
        
        Cmat = zeros(3, 3); Smat = Cmat;
        Cmat(1, 1) = 1; Smat(1, 1) = 0;
        n_max = 0; normalized = 1; Re = 1;

        [~, ~, ddU1] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r1, Re, 1-mu, ...
                                                    normalized);

        [~, ~, ddU2] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r2, Re, mu, ...
                                                    normalized);
        T = ddU1 + ddU2;
        y0 = [T(1,1); T(1, 2); T(1, 3); T(2, 2); T(2, 3); T(3, 3)];
        y = y0 + noise;
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