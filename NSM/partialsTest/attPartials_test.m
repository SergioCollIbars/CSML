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

Iner  = diag([3, 2, 1]);    % Inertia matrix;
measMode = "1";

% % % attitude error. Linear
% % Amp     = 1E-11;                                     % [rad]
% % At      = Amp.*t.^2.*ones(3, Nt).*[1;0.5;0.7];          % [rad] [yaw, pitch, roll]
% % dA_dt   = 2*Amp*t.*ones(3, Nt).*[1;0.5;0.7];
% % ddA_ddt = 2*Amp*ones(3, Nt).*[1;0.5;0.7];  

% attitude error. random
Amp     = 1E-11;                                     % [rad]
At      = normrnd(0, Amp, [3, Nt]);       % [rad] [yaw, pitch, roll]
dA_dt   = normrnd(0, Amp, [3, Nt]);
At    = Amp.*ones(3,Nt);
dA_dt = zeros(3, Nt);

% attitude nominal value
Amp = 1E-5;                                            % [rad]
attitude  = Amp.*sin(t).*ones(3, Nt).*[1;1;1];           % nominal attitude [rad]
datt_dt   = Amp.*t.*cos(t).*ones(3, Nt).*[1;1;1];

att_Err = [At;dA_dt];
[angVel_true, angAcc_true] = compute_angularVals(attitude, datt_dt, Iner, measMode, att_Err);
[angVel_nom, angAcc_nom]   = compute_angularVals(attitude, datt_dt, Iner, measMode, zeros(6, Nt));
angVel_error = angVel_true - angVel_nom;


% compute approximation error
Err_angVel = ones(3, Nt) * NaN;  Err_angAcc = ones(3, Nt) * NaN; 
for j = 3:Nt-2
    % first order partials
   if(measMode == "2")
       % GG + STT
       [~, ~, H_omega_dA, H_omegaDot_dA, H_omega_dAdt, H_omegaDot_dAdt] = ...
           compute_angularDyadPartials(angVel_nom(:, j), ...
            attitude(:, j), datt_dt(:, j),zeros(3, 1), Iner);

        % error angular velocity
        dY = angVel_true(:, j) - angVel_nom(:, j);
        Err_angVel(:, j) = dY - [H_omega_dA,H_omega_dAdt] * [At(:, j); dA_dt(:, j)];

        % error angular acceleration
        dY = angAcc_true(:, j) - angAcc_nom(:, j);
        Err_angAcc(:, j) = dY - [H_omegaDot_dA,H_omegaDot_dAdt] * [At(:, j); dA_dt(:, j)];
   elseif(measMode == "1")
       % GG + STT + GYRO
       [Hrot_omega_dyad, H_omegaDot_dyad, H_omega, H_omegaDot] = ...
           compute_angularDyadPartials_v2(angVel_nom(:, j), Iner);

       % error angular velocity
        dY = angVel_true(:, j) - angVel_nom(:, j);
        Err_angVel(:, j) = dY - H_omega * angVel_error(:,j);

        % error angular acceleration
        dY = angAcc_true(:, j) - angAcc_nom(:, j);
        Err_angAcc(:, j) = dY - H_omegaDot * angVel_error(:,j);
   end    
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



function [df] = centerDiff(f, dt, Nt)
    df = f * NaN;
    for j = 2:Nt-1
        df(:, j) = (f(:, j+1) - f(:, j-1))./(2*dt);
    end
end
