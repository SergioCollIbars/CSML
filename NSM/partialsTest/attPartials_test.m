clear;
clc;
close all;

%% TEST ATTIUDE PARTIALS


% time vector
f = 1/1;
tmax = 100;
t = linspace(0, tmax, tmax * f);
dt = t(2) - t(1);
Nt = length(t);

% attitude error. Linear
Amp     = 1E-11;                                     % [rad]
At      = Amp.*t.*ones(3, Nt).*[1;0.5;0.7];         % [rad] [yaw, pitch, roll]
dA_dt   = Amp*ones(3, Nt).*[1;0.5;0.7];
ddA_ddt = zeros(3, Nt);   

% attitude nominal value
Amp = 0;                                            % [rad]
attitude  = Amp.*t.*ones(3, Nt).*[0;1;0];           % nominal attitude [rad]
datt_dt   = Amp.*ones(3, Nt).*[0;1;0];
ddatt_ddt = zeros(3, Nt); 

[angVel_true, angAcc_true] = compute_angularVals(attitude + At, datt_dt + dA_dt, ddatt_ddt + ddA_ddt);
[angVel_nom, angAcc_nom]   = compute_angularVals(attitude, datt_dt, ddatt_ddt);

% compute approximation error
for j = 1:Nt
    % first order partials
   [~, ~, H_omega_dA, H_omegaDot_dA, H_omega_dAdt, H_omegaDot_dAdt] = ...
       compute_angularDyadPartials(flipud(angVel_nom(:, j)), ...
        attitude(:, j), datt_dt(:, j), ddatt_ddt(:, j));

    % error angular velocity
    dY = angVel_true(:, j) - angVel_nom(:, j);
    Err_angVel(:, j) = dY - [H_omega_dA,H_omega_dAdt] * [At(:, j); dA_dt(:, j)];


    % error angular acceleration
    dY = angAcc_true(:, j) - angAcc_nom(:, j);
    Err_angAcc(:, j) = dY - [H_omegaDot_dA,H_omegaDot_dAdt] * [At(:, j); dA_dt(:, j)];
end

figure()
subplot(1, 2, 1)
plot(t, Err_angVel, 'LineWidth', 2)
title('Angular velocity Taylor error')

subplot(1, 2, 2)
plot(t, Err_angAcc, 'LineWidth', 2)
title('Angular acceleration Taylor error')


figure()
subplot(1, 2, 1)
plot(t, angVel_true - angVel_nom, 'LineWidth', 2)
title('Angular velocity difference')

subplot(1, 2, 2)
plot(t, angAcc_true - angAcc_nom, 'LineWidth', 2)
title('Angular acceleration difference')
