clear;
clc;
close all;

%% INSTANT MEASUREMENT INFORMATION LIMIT.
% Description: Compute the IMIL for a point mass at different altitudes


% planet gravity constants
GM  = 4.902800118E12;            % [m^3 s^-2]
Ref = 1737.4E3;                  % [m]

% orbit altitudes
N   = 5000;
h   = linspace(10E3, 3500E3, N);  

% Measurement accuracy
sigma_m = 1E-12;                % [s^-2]

maxP = ones(1, N) * NaN;
for k = 1:N
    r    = Ref + h(k); % [m]
    maxP(k) =  r^4 * sigma_m / (3 * GM);
end

figure()
semilogx(h./1E3, maxP, 'LineWidth', 2, 'Color','r')
grid on;
xlabel('Altitude [Km]')
ylabel('[m]')
title('Instant Measurement Limit. Moon altitude')
