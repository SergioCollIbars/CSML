clc;
close all;
format long g;
clear;

%%     TEST TAYLOR SERIES EXPANSION TO SIGNAL STRENGTH

% IMPORTS
addpath('functions/');
set(0,'defaultAxesFontSize',16);

G = 6.67430e-11; % [N m^2 Kg^-2]

% analytical function. Taylor series
syms r dr
R = r + dr;
f = 1 / (R*R);

T = taylor(f, dr, 'Order', 5);

%%
% Bennu parameters
GM_Bennu = 5.2; % [N m^2 Kg^-1]
Rmin_Bennu = 250; % [m]
a_Bennu = 1000; % [m]

% Eros parameters
GM_Eros = 4.5960443E5; % [N m^2 Kg^-1]
Rmin_Eros = 16000; % [m]
a_Eros = Rmin_Eros + 40e3;   % [m]

% Earth parameters 
GM_Earth = 3.986004418e14; % [N m^2 Kg^-1]
Rmin_Earth = 6378e3; % [m]
a_Earth = Rmin_Earth + 229e3;   % [m]

% Computing Taylor series order vs the Orbit radius
nmax = 7;
GM = GM_Bennu;
Rmin = Rmin_Bennu;
points = 1000;
R = linspace(Rmin, 20 * Rmin, points);

dUa = zeros(nmax, points);
dUb = zeros(nmax, points);

Ara = 1;
Arb = -1;
for n = 1:nmax
    if (rem(n, 2) == 0)
        dUa(n, :) = NaN;
        dUb(n, :) = NaN;
    else
        for p = 1:points
            Rr = R(p);
            dUa(n, p) = (-1)^n * (n+1) * GM / (Rr^(n+2)) * Ara^n;
            dUb(n, p) = (-1)^n * (n+1) * GM / (Rr^(n+2)) * Arb^n;
        end
    end
end
ind = find(~isnan(dUa(:, 1)));
acc_diff = abs(dUa(ind, :) - dUb(ind, :));
acc_err = abs(acc_diff(1, :) - sum(acc_diff))./ sum(acc_diff);

tt = [];
for k = 1:nmax
    if(rem(k, 2) ~= 0)
        tt = [tt, "Order = " + string(k + 1)];
    end
end
tt = [tt, "Threshold"];

% compute 

% Plot 
body = "Bennu";
figure()
semilogy(R, acc_diff, 'LineWidth', 2)
xlabel('Orbital radius [m]')
ylabel(' | \partial U_A / \partial r - \partial U_B / \partial r |')
title('Gravity signal strength for different Taylor series orders. ' + string(body))
hold on;
semilogy(R, ones(1, points) * 1e-15, 'LineWidth', 2, 'LineStyle','--')
legend(tt)

figure()
plot(R, acc_err * 100, 'LineWidth', 2, 'Color', 'g')
xlabel('Orbital radius [m]')
ylabel(' Relative error [%]')
title('Gravity signal strength relative error. ' + string(body))

%%
% Compute signal resolution
Tprep = 1;  % [sec]
T = linspace(0.01, 100, points);   % [sec]

D = zeros(3, points);
DGOCE = D;
aa = [a_Earth; a_Eros; a_Bennu];
Rr = [Rmin_Earth; Rmin_Eros; Rmin_Bennu];
GMm = [GM_Earth; GM_Eros; GM_Bennu];

for i = 1:3
    a = aa(i);
    R = Rr(i);
    GM = GMm(i);
    f_orbit = sqrt(GM / (a^3)) * 1/(2*pi);
    D(i, :) = 2 * pi * R * (Tprep + 2*T) * f_orbit / 1000;
    DGOCE(i, :) = 2 * pi * R * f_orbit / 0.1 / 1000;
end

Dg_Earth = 2 * pi * Rmin_Earth / 1000 / 700;    % [Km]
Dg_Bennu = 2 * pi * Rmin_Bennu / 1000 / 9;      % [Km]
Dg_Eros = 2 * pi * Rmin_Eros / 1000 / 16;       % [Km]

% Plot
figure()
loglog(T, D, 'LineWidth', 2)
xlabel('Interrogation time [s]')
ylabel('Spatial resolution D [km]')
title('Spatial resolution for QGG vs Interrogation time')

hold on;
loglog(T, DGOCE(1, :), 'LineWidth', 2, 'LineStyle', '-.', ...
    'Color', "#0072BD")
loglog(T, DGOCE(2, :), 'LineWidth', 2, 'LineStyle', '-.', ...
    'Color', "#D95319")
loglog(T, DGOCE(3, :), 'LineWidth', 2, 'LineStyle', '-.', ...
    'Color', "#EDB120")

hold on;
loglog(T, ones(1, points) * Dg_Earth, 'LineStyle',':', 'Color', "#0072BD" ,...
    "LineWidth", 2)
loglog(T, ones(1, points) * Dg_Eros, 'LineStyle',':', 'Color', "#D95319", ...
    "LineWidth", 2)
loglog(T, ones(1, points) * Dg_Bennu, 'LineStyle',':', 'Color', "#EDB120", ...
    "LineWidth", 2)

legend('QGG Earth', 'QGG Eros', 'QGG Bennu', ...
    'GOCE Earth', 'GOCE Eros', 'GOCE Bennu', ...
    'Earth g res n = 100', 'Eros g res n = 16' , 'Bennu g res n = 9', ...
    'Location', 'northwestoutside');


