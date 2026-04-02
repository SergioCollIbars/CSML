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
normalized = file(3);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);

%% MOON VALUES FROM GRMGM 12000 FILE
R_M      = 1.7380000000000000e+06;   % [m]
GM_moon  = 4.9028001224453001e+12;   % [m^3/ s^2]
GM_sigma = 6.4536052689015518e-24;   % [m^3/ s^2]

%% GG measurements (true + perturbed)
Monte_Carlo = 400;  % number of monte carlo realizations
NH          = 2;   % number of altitudes
NS          = 4;  % number of points in the sphere
n_max       = 10; [Nc, Ns, Ncs]    = count_num_coeff(1200);
altitudes   = linspace(1E3, 100E3, NH);


% atittude error
att_err = 0.5 * 1E-3 * pi /( 3600 * 180); % [radians ]

mu_h_xx = nan(NH, 2); mu_h_yy = nan(NH, 2); mu_h_zz = nan(NH, 2);
sigma_h_xx = nan(NH, 2); sigma_h_yy = nan(NH, 2); sigma_h_zz = nan(NH, 2);
for h = 1:NH
    disp('Computing altitude = ' + string(h) + '/' + string(NH));
    Y_true      = nan(9, NS, Monte_Carlo);
    Y_nom       = nan(9, NS, Monte_Carlo);
    Y_nom_ECEF  = nan(3, 3, NS);

    disturbance_3 = nan(3, NS, Monte_Carlo);

    % Moon sphere
    [x,y,z] = fibonacci_sphere(NS, R_M + altitudes(h));

    % compute nominal measurements
    [Cnm, Snm] = list2mat(1200, Nc, Ns, SH_coeff);
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
                                    noise_att(2), 0, [3, 2, 1]); % [BN]
    
        [Cnm, Snm] = list2mat(1200, Nc, Ns, SH);
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
    D_XX  = squeeze(disturbance_2(1, :, :));          % NS x MC
    D_YY  = squeeze(disturbance_2(2, :, :));          % NS x MC
    D_ZZ  = squeeze(disturbance_2(3, :, :));          % NS x MC

    for k = 1:NS
        histogram(D_XX(k, :)); hold on;
    end
    
    D3_XX = squeeze(disturbance_3(1, :, :));
    D3_YY = squeeze(disturbance_3(2, :, :));
    D3_ZZ = squeeze(disturbance_3(3, :, :));
    
    M_D_XX = nan(1, NS); S_D_XX = nan(1, NS);
    M_D_YY = nan(1, NS); S_D_YY = nan(1, NS);
    M_D_ZZ = nan(1, NS); S_D_ZZ = nan(1, NS);

    M_D3_XX = nan(1, NS); S_D3_XX = nan(1, NS);
    M_D3_YY = nan(1, NS); S_D3_YY = nan(1, NS);
    M_D3_ZZ = nan(1, NS); S_D3_ZZ = nan(1, NS);
    
    % mean and std for MC distribution per point.
    for i = 1:NS
        M_D_XX(i) = mean(D_XX(i, :));
        M_D_YY(i) = mean(D_YY(i, :));
        M_D_ZZ(i) = mean(D_ZZ(i, :));

        S_D_XX(i) = std(D_XX(i, :));
        S_D_YY(i) = std(D_YY(i, :));
        S_D_ZZ(i) = std(D_ZZ(i, :));

        M_D3_XX(i) = mean(D3_XX(i, :));
        M_D3_YY(i) = mean(D3_YY(i, :));
        M_D3_ZZ(i) = mean(D3_ZZ(i, :));

        S_D3_XX(i) = std(D3_XX(i, :));
        S_D3_YY(i) = std(D3_YY(i, :));
        S_D3_ZZ(i) = std(D3_ZZ(i, :));
    end

    % overall mean (equally weight)
    MO_D_XX = mean(M_D_XX);
    MO_D_YY = mean(M_D_YY);
    MO_D_ZZ = mean(M_D_ZZ);

    MO_D3_XX = mean(M_D3_XX);
    MO_D3_YY = mean(M_D3_YY);
    MO_D3_ZZ = mean(M_D3_ZZ);

    % Total total variance law
    SO_D_XX = sqrt( mean(S_D_XX.^2) + mean((M_D_XX - MO_D_XX).^2) );
    SO_D_YY = sqrt( mean(S_D_YY.^2) + mean((M_D_YY - MO_D_YY).^2) );
    SO_D_ZZ = sqrt( mean(S_D_ZZ.^2) + mean((M_D_ZZ - MO_D_ZZ).^2) );

    SO_D3_XX = sqrt( mean(S_D3_XX.^2) + mean((M_D3_XX - MO_D3_XX).^2) );
    SO_D3_YY = sqrt( mean(S_D3_YY.^2) + mean((M_D3_YY - MO_D3_YY).^2) );
    SO_D3_ZZ = sqrt( mean(S_D3_ZZ.^2) + mean((M_D3_ZZ - MO_D3_ZZ).^2) );

    mu_h_xx(h, 1) = MO_D_XX;
    mu_h_yy(h, 1) = MO_D_YY;
    mu_h_zz(h, 1) = MO_D_ZZ;

    sigma_h_xx(h, 1) = SO_D_XX;
    sigma_h_yy(h, 1) = SO_D_YY;
    sigma_h_zz(h, 1) = SO_D_ZZ;


    mu_h_xx(h, 2) = MO_D3_XX;
    mu_h_yy(h, 2) = MO_D3_YY;
    mu_h_zz(h, 2) = MO_D3_ZZ;

    sigma_h_xx(h, 2) = SO_D3_XX;
    sigma_h_yy(h, 2) = SO_D3_YY;
    sigma_h_zz(h, 2) = SO_D3_ZZ;
end

% plot SH XX, YY, ZZ disturbances
plot_disturbance(mu_h_xx(:, 1), sigma_h_xx(:, 1), altitudes, "b", "SH XX component")
plot_disturbance(mu_h_yy(:, 1), sigma_h_yy(:, 1), altitudes, "b", "SH YY component")
plot_disturbance(mu_h_zz(:, 1), sigma_h_zz(:, 1), altitudes, "b", "SH ZZ component")

plot_disturbance(mu_h_xx(:, 2), sigma_h_xx(:, 2), altitudes, "g", "FRM XX component")
plot_disturbance(mu_h_yy(:, 2), sigma_h_yy(:, 2), altitudes, "g", "FRM YY component")
plot_disturbance(mu_h_zz(:, 2), sigma_h_zz(:, 2), altitudes, "g", "FRM ZZ component")


%% FUNCTIONS
function [] = plot_disturbance(mu_h, sigma_h, h, cl, tt)
    figure(); hold on; grid on;
    
    x = h./1e3;
% %     mu  = mu_h(:)';    
    sig = sigma_h(:)';
    
% %     ylow  = mu - sig;
% %     yhigh = mu + sig;
    
% %     x_fill = [x, fliplr(x)];
% %     y_fill = [ylow, fliplr(yhigh)];
    
    set(gca,'YScale','log')
    
% %     fill(x_fill, y_fill, cl, ...
% %         'FaceAlpha', 0.20, ...
% %         'EdgeColor', 'none'); hold on;
    
    plot(x, sig, 'Color', cl, 'LineWidth', 1.2);

    xlabel('Altitude [km]');
    ylabel('Disturbance');
    title(tt)
    
    grid on; xlabel('Altitude [Km]'); ylabel('milli-Eotvos');
    
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

end