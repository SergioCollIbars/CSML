clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);
addpath("functions/")
%%          STM COMPUTATION ERROR ANALYSIS
% Description: Using the gradiometer measurements, analyize the STM error
% under different dynamical systems. The STM is compared against the object
% obtained with non-linear integrator.
% Date: 04/22/2024


% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 10*86400;     % [s]
frec = 10/10;        % [Hz]
N = tmax*frec;
TIME = linspace(tmin, tmax, N);

intOrder = 4;
system = "CR3BP_inertial";     % 2BP / CR3BP / CR3Bp_inertial

% initial conditions
if(system == "2BP")
    R =  6.3781E6;  % [m]
    m = 5.9722E24;  % [Kg]
    GM = G*m;       % [m^3 s^-2]
    n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
    nmax = 3;
    planetParams = [GM, R, nmax, 0];
    poleParams = [2*pi/86400, deg2rad(180+45),...
        deg2rad(164.11), deg2rad(89.374)];
    Cmat = [1, 0, 0, 0;...
         0, 0, 0, 0;...
        -1.08E-3, 0, 1.57E-6, 0;...
         2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 
    
    Smat = [0, 0, 0, 0;...
         0, 0, 0, 0;...
         0, 0, -9.03E-7, 0;...
         0, 2.68E-7, -2.12E-7, 1.98E-7]; 

    e = 0.2;                     % eccentrycity
    a = R + 1000e3;              % semi major axis [m]
    rho = a * (1 - e^2);         % orbital param [m]
    
    i = deg2rad(45);             % inclination [rad]
    omega = deg2rad(0);          % arg periapsis [rad]
    Omega = deg2rad(0);          % RAAN [rad]
    f = deg2rad(-180);            % true anomaly [rad]
    
    [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
elseif(system == "CR3BP" || system == "CR3BP_inertial")
    R = 384399e3;               % [m]
    m_1 = 5.974E24;             % [Kg]
    m_2 = 7.348E22;             % [Kg]
    GM = G*(m_1 + m_2);         % [m^3 s^-2]
    mu =  m_2 / (m_1 + m_2);    % mass ratio
    planetParams = [mu, R, 1, 1];
    poleParams = zeros(6, 1);
    n = sqrt(GM / R^3);         % mean motion circular orbit [1/s]
    TIME = TIME.*n;             % time non-dimensionalization [-]
    planetParams(3) = n;
    Cmat = [];
    Smat = [];

    r0 = [1.021968177072928; 0; -0.18206]; 
    v0 = [0; -0.1031401430288178; 0]; % L1 orbit
end
X0 = [r0;v0];
STM_0 = reshape(eye(6,6), [36, 1]);

if(system == "CR3BP_inertial")
    [r0, v0] = rotate2inertial(X0(1:3), X0(4:6), 0, 1);
    X0 = [r0;v0];
end

% integrator
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, system), TIME, [X0; STM_0], options);

% compute gradiometer measurements
T = ones(9, N) * NaN;
if(system == "2BP")
    R_ACI = state(:, 1:3);
    [~, ddU_ACI] = gradiometer_meas(TIME, R_ACI, Cmat, Smat, ...
        poleParams, planetParams);
    T = reshape(ddU_ACI, [9, N]);
elseif(system == "CR3BP")
    T = ones(9, N) * NaN;
    for j = 1:N
        x = state(j, 1:6);
        r1 = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
        r2 = sqrt((x(1) + mu - 1)^2 + x(2)^2 + x(3)^2);
        ddU = gradmeas_rotFrame(mu, x(1), x(2), x(3),r1, r2);
        T(:, j) = reshape(ddU, [9, 1]);
    end
elseif(system == "CR3BP_inertial")
    T = ones(9, N) * NaN;
    for j = 1:N
        x = state(j, 1:6);
        ddU = gradmeas_CR3BP_inertial(mu,x(1), x(2),...
                    x(3), t(j));
        T(:, j) = reshape(ddU, [9, 1]);
    end
end

% compute reconstructed STM
[~, STM_est] = reconstruct_traj(TIME, T, X0,...
    system, intOrder);

% compute froebius norm
FN_true = ones(N, 1) * NaN;
FN_est  = ones(N, 1) * NaN;
FN_diff  = ones(N, 1) * NaN;
errC = state(:, 7:end) - STM_est;
for j = 1:N
    c = reshape(state(j, 7:end), [6,6]);
    ct = ctranspose(c);
    FN_true(j) = sqrt(trace(c* ct));

    c = reshape(STM_est(j, :), [6,6]);
    ct = ctranspose(c);
    FN_est(j) = sqrt(trace(c* ct));


    c = reshape(state(j, 7:end), [6,6]) - reshape(STM_est(j, :), [6,6]);
    ct = ctranspose(c);
    FN_diff(j) = sqrt(trace(c* ct));
end

% plot froebius norm for STM
lw = 3;
figure()
subplot(3, 1, 1)
if(system == "2BP")
    TIME = TIME./86400;
    xlg = "t [days]";
else
    xlg = "t [-]";
end
semilogy(TIME, FN_true, 'LineWidth', lw, 'Color', 'g')
xlabel(xlg)
ylabel('|STM true|')

subplot(3, 1, 2)
semilogy(TIME, FN_est, 'LineWidth', lw, 'Color', 'g')
xlabel(xlg)
ylabel('|STM est|')

subplot(3, 1, 3)
semilogy(TIME, FN_diff, 'LineWidth', lw, 'Color', 'r')
xlabel(xlg)
ylabel('|STM true - STM est|')

sgtitle('STM froebius norm')

% plot STM components and error
figure()
for j = 1:36
    subplot(6, 6, j)
    semilogy(TIME, abs(state(:, 6+j)), 'LineWidth', lw, 'Color', 'r')
    xlabel(xlg)
    ylabel('STM_{' + string(j) + '}');
end
sgtitle('True STM time series')

figure()
for j = 1:36
    subplot(6, 6, j)
    semilogy(TIME, abs(state(:, 6+j) - STM_est(:, j)), 'LineWidth', lw, ...
        'Color', 'k')
    xlabel(xlg)
    ylabel('STM_{' + string(j) + '}');
end

sgtitle('STM error time series')