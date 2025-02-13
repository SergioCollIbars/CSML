clear;
clc;
close all;

%%              STM INTEGRATION METHODS
% Description: analyze different STM integration mehods at different
% frequnecy ranges. 

set(0,'defaultAxesFontSize',16);
addpath("functions/")
addpath("../../QGG/data_files/")

% Simulation params
tmin = 0;           % [s]
tmax = 6*86400;   % [s]
frec = [0.1, 0.5, 1, 2, 5];
Nfrec = length(frec);
intOrder = [0, 2, 3, 4];

% Planet params
G = 6.67430e-11;    % [N m^2 Kg^-2]
normalized = 1;
environment = "CR3BP";     % 2BP / CR3BP

R = 384399e3;  % [m]
m_1 = 5.974E24;  % [Kg]
m_2 = 7.348E22;  % [Kg]
GM = G*(m_1 + m_2); % [m^3 s^-2]
mu = m_2 / (m_1 + m_2); % mass ratio

n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
planetParams = [mu, R, n, normalized];
poleParams = zeros(5, 1);

% initial conditions
% % X0 = [1.021968177072928; 0; -0.18206; 0; -0.1031401430288178; 0]; % L1 orbit
X0 = [1.0953533743235189E+0; -1.0879975950267760E-28; 0;...
     1.3016066486537214E-15; 2.9531900698678965E-1; 0]; % L2 orbit

% output data
error_pos = ones(length(intOrder), Nfrec, 3);
error_vel = ones(length(intOrder), Nfrec, 3);

% loop
for j = 1:length(intOrder)
    disp('computing order = ' + string(intOrder(j)))
    for i = 1:Nfrec
        % time definition
        N = tmax*frec(i);
        t = linspace(tmin, tmax, N);
        t = t.*n;           % time non-dimensionalization [-]

        % compute states and errors
        [state_true, state_recons, state_error] = orbit_corrector...
        (X0, planetParams, poleParams ,[], [], t, environment, intOrder(j));
        
        state_error  = abs(state_error);

        error_pos(j, i, 1) = max(state_error(:, 1));
        error_pos(j, i, 2) = max(state_error(:, 2));
        error_pos(j, i, 3) = max(state_error(:, 3));

        error_vel(j, i, 1) = max(state_error(:, 4));
        error_vel(j, i, 2) = max(state_error(:, 5));
        error_vel(j, i, 3) = max(state_error(:, 6));
    end
end

% plot error in position
figure()
for j = 1:3
    subplot(1, 3, j)
    variable = zeros(length(intOrder), Nfrec);
    variable = error_pos(:, :, j);
    loglog(frec, variable', 'LineWidth', 2.5, 'Marker', 'square')
    xlabel('Frequency [Hz]')
    ylabel('max error [-]')
    title('R' + string(j))
    if(j == 1)
        legend('order 0', 'order 2', 'order 3', 'order 4');
    end
    grid on;
end
sgtitle('L1 position maximum error vs frequency')


% plot error in velocity
figure()
for j = 1:3
   subplot(1, 3, j)
    variable = zeros(length(intOrder), Nfrec);
    variable = error_vel(:, :, j);
    loglog(frec, variable', 'LineWidth', 2.5, 'Marker', 'square')
    xlabel('Frequency [Hz]')
    ylabel('max error [-]')
    title('V' + string(j))
    if(j == 1)
        legend('order 0', 'order 2', 'order 3', 'order 4');
    end
    grid on;
end
sgtitle('L1 velocity maximum error vs frequency')