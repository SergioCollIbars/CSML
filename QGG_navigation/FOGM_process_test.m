clear;
clc;
close all;

f = 1/30;
tmax = 600 * 86400;
N = round(tmax * f);
TIME = linspace(0, tmax, N);
dt = TIME(2) - TIME(1);
sigmaQ = 1E-8;                                       % [E / sqrt(sec)]
varQ   = sigmaQ^2;                                   % [E^2 / sec]
q      = varQ * 1E-18;                               % [1 / sec^5]
tau    = 1*86400;                                    % [sec]

% % % compute frequency cut
% % S =(sqrt(2) * 1E-12)^2;                              % [E^2 / Hz]
% % fc= 1E-3;
% % q = S * (1+(2*pi*fc*tau)^2) / (tau^2);               % [E^2 / s]

% bais and covariance
s = q * tau / 2 * (1 - exp(-2/tau * dt));
b = zeros(1, N); 

noise = normrnd(0, 1E-12 * sqrt(f), [1, N]);
for j = 1:N-1
    w = normrnd(0, sqrt(s), [1, 1]);
    b(j + 1)    = exp(-dt/tau) * b(j) + w;      % [1 / s^2]
end

figure()
plot(TIME./3600, b./1E-9 + noise./1E-9, 'LineWidth', 2)
ylabel('Eotvos')
xlabel('TIME [hr]')
title('FOGMP bias + noise time domain')
legend('FOMP \tau = ' + string(tau));

% compute PSD of the signal residual
res = b + noise;
nfft     = [];                      % [] lets pwelch choose (or specify e.g. 2048)
window   = hamming(round(length(res)/8)); % or other window
noverlap = round(length(window)/2); % 50% overlap typical

% --- PSD estimation ---
[PSD, f] = pwelch(res, window, noverlap, nfft, f);

% --- Plot ---
figure;
loglog(f, PSD);
grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [(unit^2/Hz)]'); 
title('Power Spectral Density of Residual');

% % % integrate the STM
% % options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% % [~, STATE] = ode113(@(t, x) integrate_STM(t, x, tau), TIME, reshape(eye(2,2), [4, 1]), options);
% %  
% % PHI11 = STATE(:, 1)';
% % PHI21 = STATE(:, 2)';
% % PHI12 = STATE(:, 3)';
% % PHI22 = STATE(:, 4)';


%%  FUNCTIONS
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

function [dx] = integrate_STM(t, x, tau)
    S = reshape(x, [2,2]);
    J = [0, 1;0, 0];
    PHIdot = J * S;
    dx = reshape(PHIdot, [4, 1]);
end