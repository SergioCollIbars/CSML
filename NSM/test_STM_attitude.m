clear;
clc;
close all;



TIME   = linspace(0, 100, 1E3);
state0 = [1E-3; 1E-3; 1E-3; 1E-2; 1E-4; 1E-3];
Iner   = [2, 1, 0;1, 1, 0; 1, 0, 3];
delta0 = state0 * 1E-3; 

% integrate motion
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) att_integrator(t, x, Iner), TIME, [state0; PHI0], options);
[~, state_n] = ode113(@(t, x) att_integrator(t, x, Iner), TIME, [state0 + delta0; PHI0], options);

STM = state_t(:, 7:end);

figure()
subplot(1, 2, 1)
plot(TIME, state_t(:, 1:3), 'LineWidth', 2)
title('Euler angles')
subplot(1, 2, 2)
plot(TIME, state_t(:, 4:6), 'LineWidth', 2)
title('Angular velocity')

% approximate deviations using STM
delta_t   = state_t(:, 1:6)' - state_n(:, 1:6)';
delta_num = delta_t.*0;
for j = 1:length(TIME)
    PHI = reshape(STM(j, :), [6,6]);
    delta_num(:, j) =  PHI * delta_t(:, 1);
end

% compute numerical STM
[STM_num] = integrateSTM_attitude(TIME, state_t(:, 1:6)', 4, Iner);

% compute relative error STM
err_STM = (STM - STM_num)./STM;

% plot deviation approximation
err_delta = (delta_t - delta_num)./delta_t;

figure()
plot(TIME, err_delta(1:3, :), 'LineWidth', 2)
title('Euler angles relative error')

figure()
plot(TIME, err_delta(4:6, :), 'LineWidth', 2)
title('Angular values relative error')

%%         FUNCTIONS 
function [dx] = att_integrator(t, x, Iner)
    theta = x(2); phi = x(3);
    w = x(4:6);
    w1 = x(4); w2 = x(5); w3 = x(6);
    
    PHI =reshape(x(7:end), [6,6]);

    A  = [-sin(theta), 0, 1;...
                sin(phi)*cos(theta), cos(phi), 0;...
                cos(phi)*cos(theta), -sin(phi), 0];

    qDot = (A)\w;
    B = Iner * w;
    wDot  = (Iner)\(- cross(w, B));
    
    PsiDot   = 1/cos(theta) * (sin(phi)*w2 + cos(phi)*w3);
    ThetaDot = 1/cos(theta) * (cos(phi)*cos(theta)*w2 - sin(phi)*cos(theta)*w3);
    PhiDot   = 1/cos(theta) * (w1 + sin(phi)*sin(theta)*w2 + cos(phi)*sin(theta)*w3);

    a = sec(theta)*tan(theta)*ThetaDot*cos(theta) + ...
        1/cos(theta)*(-cos(phi)*sin(theta)*w2 + sin(phi)*sin(theta)*w3);
    b = sec(theta)*tan(theta)*PhiDot*cos(theta) + ...
        1/cos(theta)*(sin(phi)*cos(theta)*w2 + cos(phi)*cos(theta)*w3);

    C = [0, sec(theta)*tan(theta)*PsiDot*cos(theta), 1/cos(theta)*(cos(phi)*w2 - sin(phi)*w3);...
         0, a, 1/cos(theta)*(-sin(phi)*cos(theta)*w2 - cos(phi)*cos(theta)*w3);...
         0, b, 1/cos(theta)*(cos(phi)*sin(theta)*w2 - sin(phi)*sin(theta)*w3)];

    % Jacobian
    qDot_wrt_omega = inv(A);

    [t1] = compute_skewMat(w);
    [t2] = compute_skewMat(Iner * w);
    omegaDot_wrt_omega = -inv(Iner) * (t1 * Iner - t2);

    % construct Jacobian
    J = [C, qDot_wrt_omega; zeros(3, 3), omegaDot_wrt_omega];
    PHI_dot = J * PHI;

    % state derivative
    dx = [qDot;wDot;reshape(PHI_dot, [36, 1])]; 
end

function [skw] = compute_skewMat(vec)
    skw = [0, -vec(3), vec(2);...
           vec(3), 0, -vec(1);...
          -vec(2), vec(1), 0];
end