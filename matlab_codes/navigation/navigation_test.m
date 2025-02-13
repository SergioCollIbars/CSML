clear;
clc;
close all;
format long g;

%%                         NAVIGATION TEST                               %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/31/2024                                                      %
%                                                                         %
%   Description: Code to test gradiometer ability to reconstruct          %
%   trajectory based on the measurements. No nominal needed               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultAxesFontSize',16);
addpath("functions/")
addpath("../../QGG/data_files/")

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 86400;     % [s]
frec = 1/10;        % [Hz]
N = tmax*frec;
t = linspace(tmin, tmax, N);
intOrder = 4;

% Planet conditions
GM = 5.2;                % [m^3/s^2]
n_max = 9;
normalized = 1;
environment = "2BP";     % 2BP / CR3BP  / 3BP
path = "HARMCOEFS_BENNU_CD_1.txt";
[Cmat, Smat, Re] = readCoeff(path);
planetParams = [GM, Re, n_max, normalized];

W = 4.06130329511851E-4;          % [rad/s]
W0 = 0;                           % [rad/s]
RA = deg2rad(86.6388);            % [rad]
DEC = deg2rad(-65.1086);          % [rad]
poleParams = [W, W0, RA, DEC];

% Initital conditions. Orbital parameters
e = 0;                       % eccentrycity
a = 600;                    % semi major axis [m]
rho = a * (1 - e^2);         % orbital param [m]

i = deg2rad(45);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(0);              % true anomaly [rad]

if(environment == "2BP")
    [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);

    % define initial conditions. Small body
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];

elseif(environment == "CR3BP" || environment == "3BP")
    % define initial conditions. L1 orbit
    X0 = [1.021968177072928; 0; -0.18206; 0; -0.1031401430288178; 0]; % L1 orbit
    X0 = [1.0953533743235189E+0; -1.0879975950267760E-28; 0;...
         1.3016066486537214E-15; 2.9531900698678965E-1; 0]; % L2 orbit
    
    % define planet parameters
    R = 384399e3;  % [m]
    m_1 = 5.974E24;  % [Kg]
    m_2 = 7.348E22;  % [Kg]
    GM = G*(m_1 + m_2); % [m^3 s^-2]
    planetParams(1) = m_2 / (m_1 + m_2); % mass ratio
    planetParams(2) = R;   % primaries distance

    % define time dimensionalization
    n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
    t = t.*n;           % time non-dimensionalization [-]
    planetParams(3) = n;
end

% compute states and errors
[state_true, state_recons, state_error, t] = orbit_corrector...
    (X0, planetParams, poleParams ,Cmat, Smat, t, environment, intOrder);

%% PLOTS

% Plot trajectory in Inertial frame
figure();
plot3(state_true(:, 1)', state_true(:, 2)', state_true(:, 3)', ...
    'LineWidth', 2, 'Color', 'b');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Real trajectory. Inertial Frame');
axis equal;
grid on;

% Plot recons trajectory in Inertial frame
figure();
plot3(state_recons(1, :), state_recons(2, :), state_recons(3, :), ...
    'LineWidth', 2, 'Color', 'g');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Reconstructed trajectory. Inertial Frame');
axis equal;
grid on;

if(environment == "CR3BP" || environment == "3BP")
    plot_CR3BP(state_true, state_recons', t)
elseif(environment == "2BP")
    plot_2BP(state_true, state_recons, t)
end


%%      FUNCTIONS
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

function plot_2BP(state_true, state_recons, t)
    left = [1, 3, 5];
    right = [2, 4, 6];
    
    % plot reconstructed trajectory
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t./86400, state_true(:, j), 'LineWidth', 2, 'Color', 'b')
        hold on;
        plot(t./86400, state_recons(j, :), 'LineWidth', 2, 'Color', 'g')
        xlabel('TIME [days]')
        ylabel('r_' + string(j))
        legend('true', 'reconstructed')
    
        subplot(3, 2, right(j));
        plot(t./86400, state_true(:, j) - state_recons(j, :)', ...
            'LineWidth', 2, 'Color', "#FF00FF")
        xlabel('TIME [days]')
        ylabel('error [m]')
    end
    sgtitle('Reconstructed trajectory components. Inertial Frame')
    
    % plot reconstructed velocity
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t./86400, state_true(:, j+3), 'LineWidth', 2, 'Color', 'b')
        hold on;
        plot(t./86400, state_recons(j+3, :), 'LineWidth', 2, 'Color', 'g')
        xlabel('TIME [days]')
        ylabel('v_' + string(j))
        legend('true', 'reconstructed')
    
        subplot(3, 2, right(j));
        plot(t./86400, state_true(:, j+3) - state_recons(j+3, :)', ...
            'LineWidth', 2, 'Color', "#FF00FF")
        xlabel('TIME [days]')
        ylabel('error [m]')
    end
    sgtitle('Reconstructed velocity components. Inertial Frame')
    
    
    % plot showing error comparison velocity
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t./86400, abs(state_true(:, j) - state_recons(j, :)'),...
            'LineWidth', 2, 'Color', "#FF00FF")
        xlabel('TIME [days]')
        ylabel('r_' + string(j))
        legend('true - recons')
    
        subplot(3, 2, right(j));
        plot(t./86400, abs(state_true(:, j+3) - state_recons(j+3, :)'), ...
            'LineWidth', 2, 'Color', "#FF00FF")
       xlabel('TIME [days]')
       ylabel('v_' + string(j))
       legend('true - recons')
    end
    sgtitle('Absolute trajectory and velocity error components. Inertial Frame')
end

function plot_CR3BP(state_true, state_nom, t)
    left = [1, 3, 5];
    right = [2, 4, 6];
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t, state_true(:, j),...
            'LineWidth', 2, 'Color', "#FF00FF")
        hold on;
        plot(t, state_nom(:, j),...
        'LineWidth', 2, 'Color', 'g')
    
        xlabel('TIME [-]')
        ylabel('r_' + string(j))
        if(j == 1)
            legend('true', 'recons')
        end
    
        subplot(3, 2, right(j));
        plot(t, state_true(:, j+3), ...
            'LineWidth', 2, 'Color', "#FF00FF")
        hold on;
        plot(t, state_nom(:, j+3), ...
            'LineWidth', 2, 'Color', 'g')
       xlabel('TIME [-]')
       ylabel('v_' + string(j))
    end
    sgtitle('True trajectory and velocity components. Inertial Frame')
    
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t./86400, abs(state_true(:, j) - state_nom(:, j)),...
            'LineWidth', 2, 'Color', "r")
        xlabel('TIME [-]')
        ylabel('r_' + string(j))
        legend('true - recons')
    
        subplot(3, 2, right(j));
        plot(t./86400, abs(state_true(:, j+3) - state_nom(:, j+3)), ...
            'LineWidth', 2, 'Color', "r")
       xlabel('TIME [-]')
       ylabel('v_' + string(j))
       legend('true - recons')
    end
    sgtitle('Trajectory and velocity error components. Inertial Frame')

end
