clear;
clc;
close all;

%%                ORBITAL ELEMENTS IN UPM
% Description: compute how the orbital elements affect to the UPM for each
%   harmonics degree.


% Gravity model
n_max = 4;
GM = 5.1335;
R = 244.94;
W = (6.32E-13)^2;   
L = 77760;


% orbit geometry and orbit deviations
r = 800;
N = 400;

% plot eccentrycity effect on orbit
e = linspace(0, 1, N);
value_e = zeros(n_max+1, N);

for n = 0:n_max
    pol = (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
    for j =1:N
        As = 4 * pi/ L;
        sigma_t = r^(3+n) / (GM * R^n) * sqrt(W) / sqrt(pol) * sqrt(As / (4*pi));
        val = 0;
        for l = 0:fix((n+4)/2)
            A = factorial(n+4) / (factorial(l)*factorial(l)*factorial(n+4-2*l));
            val = val + A*(e(j)/2)^(2*l);
        end
        sigma = r^(3+n) * val / (GM * R^n) * sqrt(W) / sqrt(pol) * sqrt(As / (4*pi));
        sigma = sqrt(2*n + 1) * sigma;
        %value_e(n+1, j) = abs((sigma_t - sigma)./sigma_t);
        value_e(n+1, j) = sigma;
    end
end

% compute delta_GM effect
delta_GM = linspace(0, GM, N);
value_deltaGM = zeros(1, N);
for j = 1:N
    As = 4 * pi/ L;
    sigma_t = r^(3+n) / (GM * R^n) * W / sqrt(pol) * sqrt(As / (4*pi));
    sigma = r^(3+n) / ((GM + delta_GM(j)) * R^n) * W / sqrt(pol) * sqrt(As / (4*pi));
    value_deltaGM(j) = abs((sigma_t - sigma)./sigma_t); 
end

% compute reference radius sensitivity
R = linspace(246, 500, N);
G = 6.67430E-11;
b = 246;
c = 254; 
elipsoidal = 4/3 * pi * b * c;
rho = 1260;
value_R1 = zeros(n_max+1, N);
value_R2 = zeros(n_max+1, N);
for n = 0:n_max
    for j = 1:N
        spherical = 4 * pi * R(j);
        value_R1(n+1, j) = G*rho*(elipsoidal)/(R(j)^n * GM^2) + ...
            n*R(j)^(n-1)/(GM * R(j)^(2*n));
        value_R2(n+1, j) = G*rho*(spherical)/(R(j)^n * GM^2) + ...
            n*R(j)^(n-1)/(GM * R(j)^(2*n));
    end
end
Raxis = R;
R = 244.94;

% compute polar gap effect
PG = linspace(0, 1.3963, N);
value_gap = zeros(n_max + 1, N);
for n = 0:n_max
    pol = (n+1)^2 *(n+2)^2 + (n+2)^2*(n+1)*n + ...
        (n+1)*n*(n^2 -0.5) + 2*n^2*(n+1) + (n+1)^2;
    A = (2*n + 1) * factorial(2*n) / (factorial(n)^2 * 2^(2*n));
    for j = 1:N
        As = 4 * pi/ L;
        val = 0;
        delta = PG(j);
        for k = 0:n
            b = nchoosek(n,k);
            val = val  + ...
                b * (-1)^k / (2*n -2*k + 1) * (2 * cos(delta)^(2*n - 2*k + 1));
        end
        val = val * (-1)^n;
        I = 2*pi * A * val;

        sigma = r^(3+n) / (GM * R^n) * sqrt(W) / sqrt(pol) * sqrt(As / (4*pi));
        sigma_PG = r^(3+n) / (GM * R^n) * sqrt(W) / sqrt(pol) * sqrt(As / I);
        value_gap(n+1, j) = abs((sigma - sigma_PG) / sigma);
        %value_gap(n+1, j) = sigma_PG;
    end
end


%% PLOTS

% plot eccenctricity effect
figure();
semilogy(e, value_e, 'LineWidth', 1.5)
val_legend = {};
for j = 0:n_max
    val_legend{end+1} = "n = " + string(j); 
end
legend(val_legend)
xlabel('Eccentrycity, e [-]')
ylabel('uncertainty, \sigma [-]')
%ylabel('Relative error[-]')
title('Eccentricity effect on the UPM')

% plot reference radius effect
figure();
lw = 1.5;
semilogy(Raxis, value_R1, 'LineWidth',lw)
xlabel('Reference Radius, R [m]')
ylabel('Normalized sensitivity')
title('Elipsoidal volume \sigma sensitivity')
legend(val_legend)

figure();
semilogy(Raxis, value_R2, 'LineWidth',lw)
xlabel('Reference Radius, R [m]')
ylabel('Normalized sensitivity')
title('Spherical volume \sigma sensitivity')
legend(val_legend)


% plot delta GM error effect
figure();
plot(delta_GM./GM, value_deltaGM, 'LineWidth', 1.5, 'Color','r');
xlabel('\delta GM / GM [-]')
ylabel('relative error [-]')
title('Uncertainty in GM value effect in the UPM')

% plot polar gap effect
figure();
semilogy(rad2deg(PG), value_gap, 'LineWidth', 1.5)
legend(val_legend)
xlabel('Polar gap value, [deg]');
%ylabel('uncertainty, \sigma [-]');
ylabel('relative error, [-]')
title('Polar gap effect in the UPM');
    