clear; clc; close all;
format long g;

addpath('../../data/'); addpath(genpath("../../functions/"));
set(0,'defaultAxesFontSize',16);
%%          SH COEFFICIENTS ERROR BUDGET PER ALTITUDE
% Description: Compute the impact of gravity field coefficient uncertainty
% at different orbit altitudes. (error budget)
% Author: Sergio Coll-Ibars
% Date: 01/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gravity field
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
R_M = file(2)*1E3;   normalized = file(3); GM_moon = 4.9028001224453001e+03 * 1E9;

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);

%% GG measurements (true + perturbed)
Monte_Carlo = 100;  % number of monte carlo realizations
NH          = 40;   % number of altitudes
NS          = 1000;  % number of points in the sphere
n_max       = 20; [Nc, Ns, Ncs]    = count_num_coeff(n_max);
altitudes   = linspace(1E3, 100E3, NH);

mean_distubance = nan(1, NH);
std_distubance  = nan(1, NH);
for h = 1:NH
    disp('Computing altitude = ' + string(h) + '/' + string(NH));
    Y_true      = nan(9, NS, Monte_Carlo);
    Y_nom       = nan(9, NS, Monte_Carlo);

    % Moon sphere
    [x,y,z] = fibonacci_sphere(NS, R_M + altitudes(h));

    % compute nominal measurements
    [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);
    for j = 1:NS
        r = [x(j); y(j); z(j)];  % ECEF coordinates

        [~, ~, T_ECEF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                         r, R_M, GM_moon, normalized);

        ENU_ECEF = ecef2enu(r);
        T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
        Y_ENU    = [T_ENU(1,1);T_ENU(1,2);T_ENU(1,3);T_ENU(2,1);...
            T_ENU(2,2);T_ENU(2,3);T_ENU(3,1);T_ENU(3,2);T_ENU(3,3)];
        Y_nom(:, j, :) = Y_ENU.*ones(9, Monte_Carlo);
    end

    for k = 1:Monte_Carlo
        err = normrnd(0, SH_uncrt(2:end));
        SH  = SH_coeff + [0;err];
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH);
        parfor j = 1:NS
            r = [x(j); y(j); z(j)];   % ECEF coordinates                          
            [~, ~, T_ECEF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                         r, R_M, GM_moon, normalized);

            ENU_ECEF = ecef2enu(r);
            T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
            Y_ENU    = [T_ENU(1,1);T_ENU(1,2);T_ENU(1,3);T_ENU(2,1);...
                T_ENU(2,2);T_ENU(2,3);T_ENU(3,1);T_ENU(3,2);T_ENU(3,3)];
            Y_true(:, j, k) = Y_ENU;
        end
    end

    % disturbance
    disturbance_1 = (Y_true - Y_nom)./1E-12;  % [milli-Eotvos]
    disturbance_2 = [disturbance_1(1, :, :);disturbance_1(5, :, :);...
        disturbance_1(9, :, :)];
    
    % disturbance_h is 9 x NS x MC (2-D norm)
    D = squeeze(vecnorm(disturbance_2, 2, 1) );   % NS x MC

    % area weigths (quasi-uniform)
    w = ones(NS, 1)./NS;
    
    % mean MC per point on the globe. 
    mu_pt  = mean(D, 2);          % NSx1

    % mean and std over the globe
    mu_global     = sum(w .* mu_pt);
    std_mu_global = sqrt( sum(w .* (mu_pt - mu_global).^2) );

    mean_distubance(h)    = mu_global;
    std_distubance(h)     = std_mu_global;
end

%%
% plot disturbance RMS value
figure(); hold on; grid on;

x = altitudes./1e3;
mu = mean_distubance;
sig = std_distubance;

% --- Shaded ±1σ band ---
x_fill = [x, fliplr(x)];
y_fill = [mu - sig, fliplr(mu + sig)];

hfill = fill(x_fill, y_fill, 'b', ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', 'none');

% --- Mean curve on top ---
semilogy(x, mu, 'b', 'LineWidth', 2);

set(gca,'YScale','log');
xlabel('Altitude [km]');
ylabel('Disturbance');

grid on; xlabel('Altitude [Km]'); ylabel('milli-Eotvos');
title('Disturbance mean +- \sigma value');

y = 1;
yline(y,'LineWidth',1.5, 'Color', 'r');
text(mean(xlim), y, 'GOCE noise level', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize',11);

y = 0.01;
yline(y,'LineWidth',1.5, 'Color', 'g');
text(mean(xlim), y, 'QUANTUM noise level', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize',11);

y = 100;
yline(y,'LineWidth',1.5, 'Color', 'c');
text(mean(xlim), y, 'MEMS noise level', ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', ...
     'FontSize',11);
