clear;
clc;
close all;
format long g;


% time interval
rev = 1;
T = 10;
f = 1/5;
t = linspace(0, rev*T, rev*T * f);
dt = t(2) - t(1);
Nt = length(t);

% attitude error
frec    = 1E-4;                                                % [rad/s]
Amp     = 5E-7;
At      = Amp.*sin(frec.*t).*ones(3, Nt);                      % [rad] [yaw, pitch, roll]
dA_dt   = Amp * frec * cos(frec*t).*ones(3, Nt);
ddA_ddt = -Amp.*frec^2.*sin(frec.*t).*ones(3, Nt);   

% attitude nominal value
Amp = 0;
attitude  = Amp.*sin(frec.*t).*ones(3, Nt);                         % nominal attitude [rad]
datt_dt   = Amp.*frec*cos(frec.*t).*ones(3, Nt);
ddatt_ddt = Amp.*-frec^2*sin(frec.*t).*ones(3, Nt); 

% compute angular values
[angVel_true, angAcc_true, ~, ~] = compute_angularVals(attitude + At, datt_dt + dA_dt, ddatt_ddt + ddA_ddt, t);
[angVel_nom, angAcc_nom, H_angVel, H_angAcc]   = ...
    compute_angularVals(attitude, datt_dt, ddatt_ddt, t);

% Angular velocity partials
[H_dA, H_dA_dT] = compute_angVel_partials(attitude, datt_dt);

angVel_diff = angVel_true - angVel_nom;
figure()
plot(t, angVel_diff)
ylabel('[rad/s]')

% Approximation test
Err1  =zeros(3, Nt); Err2 = Err1;
for j = 1:Nt
    finish = 3*j; init = finish - 2;
    h_dA = H_dA(init:finish, :);
    h_dA_dT = H_dA_dT(init:finish, :);

    Err1(:, j) = h_dA * At(:, j);
    Err2(:, j) = h_dA_dT * dA_dt(:, j) + h_dA * At(:, j);
end
figure()
plot(t, vecnorm(angVel_diff - Err1))
hold on;
plot(t, vecnorm(angVel_diff - Err2))
legend('Err1', 'Err2')

%% FUNCTIONS
function [H_dA, H_dA_dT] = compute_angVel_partials(attitude, dA_dt)
    Nt  = length(attitude(1, :));
    for j = 1:Nt                            % attitude rate of change (vel)
        theta = attitude(2, j); phi = attitude(3, j);
        A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
        
        % partials of angular velocity w.r.t Euler angles
        finish = 3*j;
        init = finish - 2;
        psiDot = dA_dt(1, j); thetaDot = dA_dt(2, j);
        B = [0, -cos(theta)*psiDot, 0;...
             0, -sin(phi)*sin(theta)*psiDot, cos(phi)*cos(theta)*psiDot-sin(phi)*thetaDot;...
             0, -cos(phi)*sin(theta)*psiDot, -sin(phi)*cos(theta)*psiDot-cos(phi)*thetaDot];
        H_dA(init:finish, :) = B; 
        H_dA_dT(init:finish, :) = A;
    end
end