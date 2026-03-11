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
n_max       = 200; [Nc, Ns, Ncs]    = count_num_coeff(n_max);
altitudes   = linspace(1E3, 100E3, NH);

% atittude error
att_err = 0.5 * 1E-3 * pi /( 3600 * 180); % [radians ]

mean_distubance = nan(NH, 2);
std_distubance  = nan(NH, 2);
for h = 1:NH
    disp('Computing altitude = ' + string(h) + '/' + string(NH));
    Y_true      = nan(9, NS, Monte_Carlo);
    Y_nom       = nan(9, NS, Monte_Carlo);
    Y_nom_ECEF  = nan(3, 3, NS);

    disturbance_3 = nan(3, NS, Monte_Carlo);

    % Moon sphere
    [x,y,z] = fibonacci_sphere(NS, R_M + altitudes(h));

    % compute nominal measurements
    [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);
    parfor j = 1:NS
        r = [x(j); y(j); z(j)];  % ECEF coordinates

        [~, ~, T_ECEF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                         r, R_M, GM_moon, normalized);

        ENU_ECEF = ecef2enu(r);
        T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
        Y_ENU    = [T_ENU(1,1);T_ENU(1,2);T_ENU(1,3);T_ENU(2,1);...
            T_ENU(2,2);T_ENU(2,3);T_ENU(3,1);T_ENU(3,2);T_ENU(3,3)];
        Y_nom(:, j, :) = Y_ENU.*ones(9, Monte_Carlo);

        Y_nom_ECEF(:, :, j) = T_ECEF;
    end

    for k = 1:Monte_Carlo
        % Gravity error
        err = normrnd(0, SH_uncrt(2:end));
        SH  = SH_coeff + [0;err];

        % attitude error
        noise_att   = normrnd(0, att_err, [2, 1]);
        R          = rotationMatrix(noise_att(1), ...
                                    noise_att(2), 0, [1, 2, 3]); % [BN]
    
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

            % attitude error
            Y_nom_mat    = Y_nom_ECEF(:, :, j);
            dist_att = (R' * Y_nom_mat * R) - Y_nom_mat;
            disturbance_3(:, j, k) = [dist_att(1, 1); dist_att(2, 2);...
                dist_att(3,3)]./1E-12; % [milli-Eotvos]
        end
    end

    % disturbance
    disturbance_1 = (Y_true - Y_nom)./1E-12;  % [milli-Eotvos]
    disturbance_2 = [disturbance_1(1, :, :);disturbance_1(5, :, :);...
        disturbance_1(9, :, :)];
    
    % disturbance_h is 9 x NS x MC (2-D norm)
    D  = squeeze(vecnorm(disturbance_2, 2, 1) );   % NS x MC
    D3 = squeeze(vecnorm(disturbance_3, 2, 1) );   % NS x MC

    % area weigths (quasi-uniform)
    w = ones(NS, 1)./NS;
    
    % mean MC per point on the globe. 
    mu_pt    = mean(D, 2);          % NSx1
    mu_pt_3  = mean(D3, 2);         % NSx1

    % mean and std over the globe
    mu_global     = sum(w .* mu_pt);
    std_mu_global = sqrt( sum(w .* (mu_pt - mu_global).^2) );

    mean_distubance(h, 1)    = mu_global;
    std_distubance(h, 1)     = std_mu_global;

    mu_global     = sum(w .* mu_pt_3);
    std_mu_global = sqrt( sum(w .* (mu_pt_3 - mu_global).^2) );

    mean_distubance(h, 2)    = mu_global;
    std_distubance(h, 2)     = std_mu_global;
end

%%
% plot disturbance RMS value
figure(); hold on; grid on;

colorPalet = ['b', 'g'];
for n = 1:2
x = altitudes./1e3;
mu  = mean_distubance(:, n)';
sig = std_distubance(:, n)';

% --- Shaded ±1σ band ---
x_fill = [x, fliplr(x)];
y_fill = [mu - sig, fliplr(mu + sig)];

hfill = fill(x_fill, y_fill, colorPalet(n), ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', 'none');

% --- Mean curve on top ---
semilogy(x, mu, colorPalet(n), 'LineWidth', 2); hold on;
end

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

% % y = 100;
% % yline(y,'LineWidth',1.5, 'Color', 'c');
% % text(mean(xlim), y, 'MEMS noise level', ...
% %      'HorizontalAlignment','center', ...
% %      'VerticalAlignment','bottom', ...
% %      'FontSize',11);
