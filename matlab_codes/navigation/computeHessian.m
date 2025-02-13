clear;
clc;
close all;

addpath('functions/')
addpath("../../QGG/data_files/")
set(0,'defaultAxesFontSize',16);
%%                  COMPUTE HESSIAN MATRIX
% Description: Look for the hessian matrix in the STM propagation

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 1*8000;   % [s]
frec = 1/10;        % [Hz]
N = tmax*frec;
t = linspace(tmin, tmax, N);
system = "2BP";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earth case
% define planet parameters
% % R =  6.3781E6;  % [m]
% % m = 5.9722E24;  % [Kg]
% % GM = G*m;       % [m^3 s^-2]
% % n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
% % nmax = 3;
% % planetParams = [GM, R, nmax, 0];
% % poleParams = [2*pi/86400, 0, -pi/2, pi/2];
% % Cmat = [1, 0, 0, 0;...
% %      0, 0, 0, 0;...
% %     -1.08E-3, 0, 1.57E-6, 0;...
% %      2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 
% % 
% % Smat = [0, 0, 0, 0;...
% %      0, 0, 0, 0;...
% %      0, 0, -9.03E-7, 0;...
% %      0, 2.68E-7, -2.12E-7, 1.98E-7]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bennu case
GM = 5.2;       % [m^3/s^2]
n_max = 9;
normalized = 1;
path = "HARMCOEFS_BENNU_CD_1.txt";
[Cmat, Smat, Re] = readCoeff(path);
planetParams = [GM, Re, n_max, normalized];

W = 4.06130329511851E-4;          % [rad/s]
W0 = 0;                           % [rad/s]
RA = deg2rad(86.6388);            % [rad]
DEC = deg2rad(-65.1086);          % [rad]
poleParams = [W, W0, RA, DEC];


if(system == "2BP")
    % Initital conditions. Orbital parameters
    e = 0.2;                     % eccentrycity
    a = 1000;             % semi major axis [m]
    rho = a * (1 - e^2);       % orbital param [m]
    
    i = deg2rad(45);            % inclination [rad]
    omega = deg2rad(0);        % arg periapsis [rad]
    Omega = deg2rad(0);        % RAAN [rad]
    f = deg2rad(-180);         % true anomaly [rad]
    
    [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
    deltaX = [100;0;0;0;0;0];
elseif(system == "CR3BP")
    R = 384399e3;  % [m]
    m_1 = 5.974E24;  % [Kg]
    m_2 = 7.348E22;  % [Kg]
    GM = G*(m_1 + m_2); % [m^3 s^-2]
    mu =  m_2 / (m_1 + m_2); % mass ratio
    planetParams = [mu, R, 1, 1];
    poleParams = zeros(6, 1);


    % define time dimensionalization
    n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
    t = t.*n;           % time non-dimensionalization [-]
    planetParams(3) = n;
    
    % define initial conditions. L1 orbit
    X0 = [1.021968177072928; 0; -0.18206; 0; -0.1031401430288178; 0]; % L1 orbi
    % % X0 = [1.0953533743235189E+0; -1.0879975950267760E-28; 0;...
    % %      1.3016066486537214E-15; 2.9531900698678965E-1; 0]; % L2 orbit
end

% define initial conditions. Small body
STM_0 = reshape(eye(6,6), [36, 1]);

% define integration options
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% ODE 113
[time, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, system), t, [X0; STM_0], options);

% plot orbit states
plotOrbit(time, state_true, state_true);

% % % look for Hessian matrix value
% % STM = state_true(:, 7:end);
% % Nt = length(time);
% % Nt = 200;
% % lambda = ones(12, Nt)*NaN;
% % determinant = zeros(1, Nt);
% % H = zeros(6*(Nt), 6*(Nt));
% % for s = 2:Nt
% %     P = reshape(STM(s, :), [6,6]);
% %     p11 = P(1:3,1:3);
% %     p12 = P(1:3,4:6);
% %     p21 = P(4:6,1:3);
% %     p22 = P(4:6,4:6);
% % 
% %     B = P'*P; 
% %     b11 = B(1:3,1:3);
% %     b12 = B(1:3,4:6);
% %     b21 = B(4:6,1:3);
% %     b22 = B(4:6,4:6);
% %     
% %     % Hessian at time s
% %     Z = [0,2,0;-2,0,0;0,0,0];
% %     if(system == "2BP")
% %         Hs = [P'*P, -P';-P, eye(6,6)];
% %     elseif(system == "CR3BP")   % WARNING. This is not correct!
% %         Hs = [b11,b12,-p11,-p21;...
% %             b21+b22*Z, b22, -p12-p22*Z, -p22;...
% %             -p11-p12*Z,-p12,eye(3,3),zeros(3,3);...
% %             -p21-p22*Z, -p22, Z, eye(3,3)];
% %     end
% % 
% %     determinant(s) = det(Hs);
% %     nd = length(Hs(:, 1));
% %     D = eig(Hs);
% %     lambda(1:nd, s) = D;
% % 
% %     % Hessian for: J = sum(E'*E)
% %     H(1:6,1:6) =  H(1:6,1:6) + P'*P;
% %     init = 6*(s-1) + 1;
% %     final = init + 5;
% %     H(1:6, init:final) = -P';
% %     H(init:final,1:6) = -P;
% %     H(init:final,init:final) = eye(6,6);
% % end

Nt = length(time);
DF = zeros(6*Nt-6, 6*Nt-6);
for j = 1:Nt-1
    tspan = [t(j), t(j+1)];
    [~, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, system), tspan, [X0; STM_0], options);

    PHI_j0 = state(end, 7:end);
    g = reshape(PHI_j0, [6,6]) - eye(6,6);
    maxId = 6*j;
    minId = maxId -5;
    DF(minId:maxId, minId:maxId) = g;

    % update initial state
    X0 = state(end, 1:6)';
end


% plot hessian matrix
figure()
for k = 1:12
    subplot(6, 2, k)
    plot(1:Nt, lambda(k, :), "LineWidth", 2, "Color", "r")
end
sgtitle("Hessian eigenvalues")

figure()
plot(1:Nt, determinant, "LineWidth", 2)
ylabel("|H|")
xlabel("time index")
title("Hessian determinant")

%%
function [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e)
    % compute initial state vectors (orbit frame)
    r0 = rho / (1 + e * cos(f)) * [cos(f);...
                                  sin(f);...
                                  0];
    v0 = sqrt(GM / rho) * [-sin(f);...
                          e + cos(f);...
                          0];
    
    % compute rotation matrix: Body to ACI
    [BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
    
    % rotate initial state vectors. {I, J, K} frame
    r0 = BN' * r0;
    v0 = BN' * v0;
end

function plotOrbit(time, state_true, stateDot_true)
    figure()
    plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
        'LineWidth', 2, 'Color', 'r');
    title('Orbit. Inertial frame')
    axis equal
    grid on;
    
    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, state_true(:, j), 'LineWidth', 2, 'Color', 'b')
        xlabel('TIME [s]')
        ylabel('r_' + string(j) + ' [m]')
    end
    sgtitle('Position. Inertial frame')
    
    
    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, state_true(:, j+3), 'LineWidth', 2, 'Color', 'g')
        xlabel('TIME [s]')
        ylabel('v_' + string(j) + ' [m/s]')
    end
    sgtitle('Velocity. Inertial frame')

    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, stateDot_true(:, j+3), 'LineWidth', 2, 'Color', ...
            "#FF00FF")
        xlabel('TIME [s]')
        ylabel('a_' + string(j) + '[m/s^2]')
    end
    sgtitle('Acceleration. Inertial frame')
end


