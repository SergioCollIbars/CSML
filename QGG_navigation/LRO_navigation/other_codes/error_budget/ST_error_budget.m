clear; clc;
close all;

set(0,'defaultAxesFontSize',16);
addpath('../../data/'); addpath(genpath("../../functions/"));
%% SOLID TIDES PERTURBATION ERROR BUDGET
% Description: Compute the maximum deviation on the n = 2 & 3 coefficients
% based on the Love numbers.
% Date: 01/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";

file             = readmatrix(input_gravField);
normalized       = file(3);

SH_coeff         = file(4:end);

%% Planet constants
GM_E   = 3.986004418e14;     % Earth GM [m^3/s^2]
GM_M   = 4.9048695e12;       % Moon GM  [m^3/s^2]
GM_S   = 1.32712440018E20;   % Sun GM  [m^2/s^2]​​
R_M    = 1738.0e3;           % Moon mean radius [m]

r_EM_min = 3.633E8;          % mis distance Moon - Earth [m]
r_SM_min = 1.471E11;         % min distance Moon - Sun [m]

%% Love Numbers from GRAIL
k20 = 0.02405;  s_k20 = 0.00018;
k21 = 0.02414;  s_k21 = 0.00025;
k22 = 0.02394;  s_k22 = 0.00028;
k30 = 0.0089;   s_k30 = 0.0021;

love_numb_list      = [k20, k21, k22, k30];
love_numb_list_unct = [s_k20, s_k21, s_k22, s_k30];

%% GG measurements (true + perturbed)
Monte_Carlo = 400;  % number of monte carlo realizations % 200
NH          = 40;   % number of altitudes   % 40
NS          = 1000;  % number of points in the sphere % 1000
n_sim       = 4; [Nc, Ns, Ncs]    = count_num_coeff(1200);
altitudes   = linspace(1E3, 100E3, NH);

mu_h_xx = nan(NH, 1); mu_h_yy = nan(NH, 1); mu_h_zz = nan(NH, 1);
sigma_h_xx = nan(NH, 1); sigma_h_yy = nan(NH, 1); sigma_h_zz = nan(NH, 1);
for h = 1:NH
    disp('Computing altitude = ' + string(h) + '/' + string(NH));
    Y_true      = nan(9, NS, Monte_Carlo);
    Y_nom       = nan(9, NS, Monte_Carlo);

    % Moon sphere
    [x,y,z] = fibonacci_sphere(NS, R_M + altitudes(h));

    % compute nominal measurements
    [Cnm_t, Snm_t] = list2mat(1200, Nc, Ns, SH_coeff);
    for j = 1:NS
        r = [x(j); y(j); z(j)];  % ECEF coordinates

        [~, ~, T_ECEF] = potentialGradient_nm(Cnm_t, Snm_t, n_sim, ...
                         r, R_M, GM_M, normalized);

        ENU_ECEF = ecef2enu(r);
        T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
        Y_ENU    = [T_ENU(1,1);T_ENU(1,2);T_ENU(1,3);T_ENU(2,1);...
            T_ENU(2,2);T_ENU(2,3);T_ENU(3,1);T_ENU(3,2);T_ENU(3,3)];
        Y_nom(:, j, :) = Y_ENU.*ones(9, Monte_Carlo);
    end

    for k = 1:Monte_Carlo
        % perturb solid tides (static)
        [love_numb] = create_love_matrix(1200, love_numb_list, ...
                    love_numb_list_unct);

        % perturb monthly variations
        [deltaC_month, deltaS_month] = Max_solid_tide_monthly_var(1200);
        
        % Max gravity field perturbation due to solid tides
        [deltaC, deltaS] = Max_solid_tidial_var(1200,n_sim, love_numb, ...
                                R_M, GM_M, GM_E, GM_S, r_EM_min, r_SM_min);

        Cnm = Cnm_t + deltaC + deltaC_month; 
        Snm = Snm_t + deltaS + deltaS_month;
        parfor j = 1:NS 
            r = [x(j); y(j); z(j)];   % ECEF coordinates                          
            [~, ~, T_ECEF] = potentialGradient_nm(Cnm, Snm, n_sim, ...
                         r, R_M, GM_M, normalized);

            ENU_ECEF = ecef2enu(r);
            T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
            Y_ENU    = [T_ENU(1,1);T_ENU(1,2);T_ENU(1,3);T_ENU(2,1);...
                T_ENU(2,2);T_ENU(2,3);T_ENU(3,1);T_ENU(3,2);T_ENU(3,3)];
            Y_true(:, j, k) = Y_ENU;
        end
    end

    disturbance_1 = (Y_true - Y_nom)./1E-12;  % [milli-Eotvos]
    disturbance_2 = [disturbance_1(1, :, :);disturbance_1(5, :, :);...
        disturbance_1(9, :, :)];

     % disturbance_h is 9 x NS x MC (2-D norm)
    D_XX  = squeeze(disturbance_2(1, :, :));          % NS x MC
    D_YY  = squeeze(disturbance_2(2, :, :));          % NS x MC
    D_ZZ  = squeeze(disturbance_2(3, :, :));          % NS x MC

    % Disturbance 2
    [M_D_XX, Q3_XX] = compute_med_Q3(D_XX);
    mu_h_xx(h, 1) = M_D_XX;
    sigma_h_xx(h, 1) = Q3_XX;

    [M_D_ZZ, Q3_ZZ] = compute_med_Q3(D_ZZ);
    mu_h_zz(h, 1) = M_D_ZZ;
    sigma_h_zz(h, 1) = Q3_ZZ;

    [M_D_YY, Q3_YY] = compute_med_Q3(D_YY);
    mu_h_yy(h, 1) = M_D_YY;
    sigma_h_yy(h, 1) = Q3_YY;
end

plot_disturbance(mu_h_xx, sigma_h_xx, altitudes, "m", "XX component")
plot_disturbance(mu_h_yy, sigma_h_yy, altitudes, "m", "YY component")
plot_disturbance(mu_h_zz, sigma_h_zz, altitudes, "m", "ZZ component")


%% FUNCTIONS
function [MEDIAN, Q3] = compute_med_Q3(data)
    % flaten data
    x = abs(data(:));

    MEDIAN = median(x);
    Q3     = prctile(x, 95);
end

function [deltaC, deltaS] = Max_solid_tidial_var(n_max,n_sim, love_numb, ...
    R_M, GM_M, GM_E, GM_S, r_E, r_S)
    
    deltaC = zeros(n_max+1, n_max+1);  deltaS = zeros(n_max+1, n_max+1);
    for n = 0:n_sim
        for m = 0:n
            row_index = n + 1; col_index = m+1;

            love_val  = love_numb(row_index, col_index);
            
            x = linspace(-1,1,20001);

            % Fully normalized associated Legendre
            Pbar = legendre(n, x, 'norm');
            Pnm  = squeeze(Pbar(m+1,:));

            [maxVal_Pnm, ~] = max(abs(Pnm));

            Earth_pert = GM_E / GM_M * (R_M / r_E)^(n+1) * maxVal_Pnm;
            Sun_pert   = GM_S / GM_M * (R_M / r_S)^(n+1) * maxVal_Pnm;

            deltaC(row_index, col_index) = love_val / (2*n + 1) * ...
                (Earth_pert + Sun_pert);
            if(m == 0)
                deltaS(row_index, col_index) = 0;
            else
                deltaS(row_index, col_index) =  love_val / (2*n + 1) * ...
                    (Earth_pert + Sun_pert);
            end
        end
    end
end

function [love_numb] = create_love_matrix(n_max, love_num_list, ...
    love_num_unct)
    love_numb = zeros(n_max+1, n_max+1);

    vals = love_num_list + love_num_unct .* randn(size(love_num_list));
     
    love_numb(3, 1:3)      = [vals(1), vals(2), vals(3)];
    love_numb(4, 1)        = vals(4);
end

function [deltaC, deltaS] = Max_solid_tide_monthly_var(n_max)
        % amplitudes
                               % C20, C21, S21, C22, S22
        amplitudes = (1E-10).*[-1.89, 0, 0, 3.28, -4.41;...
                        0, 4.72, -0.02, 0, 0;...
                    -0.36,0,0,0.58,-0.88;...
                    -0.32,0,0,0.75,-0.48;...
                    -0.16,0,0,0.51,-0.51;...
                     0,0.64,0.26,0,0];

        deltaC = zeros(n_max+1, n_max+1); deltaS = deltaC;

        deltaC(3, 1)      = sum(amplitudes(:,1)); % C20
        deltaC(3, 2)      = sum(amplitudes(:,2)); % C21
        deltaC(3, 3)      = sum(amplitudes(:,4)); % C22

        deltaS(3, 2)      = sum(amplitudes(:,3)); % S21
        deltaS(3, 3)      = sum(amplitudes(:,5)); % S22
end

function [] = plot_disturbance(mu_h, sigma_h, h, cl, tt)
    figure(); hold on; grid on;
    
    x = h./1e3;
    sig = sigma_h(:);
    
    set(gca,'YScale','log')
    
    plot(x,  sig, 'Color', cl, 'LineWidth', 1.2);

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