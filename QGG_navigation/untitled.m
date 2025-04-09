clear;
clc;
close all;

f = 1/30;
tmax = 6 * 86400;
N = round(tmax * f);
TIME = linspace(0, tmax, N);
dt = TIME(2) - TIME(1);
sigmaQ = 1E-5;                                       % [E / sec]
varQ   = sigmaQ^2;                                   % [E^2 / sec^2]
qE     = varQ*dt;                                    % [E^2 / sec] assuming it is multiplied by dt
q      = qE * 1E-18;                                 % [1 / sec^5]
tau    = 5;                                          % [sec]

% bais and covariance
s = q * tau / 2 * (1 - exp(-2/tau * dt));
b = zeros(1, N); 
b_RW = b;

% uncertainty 
P0 = eye(2,2).*(0)^2;
Sb = zeros(1, N); Sd = Sb;

noise = normrnd(0, 1E-12, [1, N]);
for j = 1:N-1
    w = normrnd(0, sqrt(s), [1, 1]);
    w2 = normrnd(0, sqrt(q*dt), [1, 1]);
    b(j + 1)    = exp(-dt/tau) * b(j) + w;      % [1 / s^2]
    b_RW(j + 1) = b_RW(j) + w2;
    

    [PHIbd] = compute_STM_FOGM(dt, tau);
    [Sbd] = compute_PN_FOGM(dt, tau, q);
    P = PHIbd * P0 * PHIbd' + Sbd;
    Sb(j) = sqrt(P(1,1)); Sd(j) = sqrt(P(2,2));

    P0 = P;
end

figure()
plot(1:N, b./1E-9 + noise./1E-9, 'LineWidth', 2, 'color', 'b')
hold all;
%plot(1:N, b_RW./1E-9 + noise./1E-9, 'LineWidth', 2)
legend('FOMP \tau = ' + string(tau), 'Random Walk')
ylabel('Eotvos')


% FUNCTIONS
function [PHIbd] = compute_STM_FOGM(t, tau)
    PHIbd = [1, tau*(1 - exp(-t/ tau));...
        0, exp(-t/tau)];
end

function [Sbd] = compute_PN_FOGM(t, tau, sigmaQ_b)
    S11 = tau^2 * ((1-exp(-2*t/tau)) - 4*(1-exp(-t/tau)) + 2*t/tau);
    S12 = tau * (1 - exp(-t/tau))^2;
    S21 = tau * (1 - exp(-t/tau))^2;
    S22 = (1 - exp(-2*t/tau));

    Sbd = sigmaQ_b * tau / 2 * [S11,S12;S21,S22];
end