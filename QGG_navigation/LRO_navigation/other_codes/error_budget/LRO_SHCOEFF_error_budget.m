clear; clc; close all;
format long g;
set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')
%%          SH COEFFICIENTS ERROR BUDGET LRO
% Description: Compute the impact of gravity field coefficient uncertainty
% in the LRO orbit. (error budget)
% Author: Sergio Coll-Ibars
% Date: 01/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRGM1200 gravity field 
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";
input_clone      = "HARMCOEFS_MOON_CLONE_01_GRGM1200PRIM.txt";
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
n_max            = file(1); normalized = file(3);
[Nc, Ns, ~]      = count_num_coeff(n_max);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_clone);
SH_clone         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);


%% MOON VALUES FROM SPICE
% % % Get GM for the Moon [m^3/s^2]
% % [GM_moon] = cspice_bodvrd('MOON', 'GM', 1)*1E9; 
% % 
% % % get Earth Radius [m]
% % [radii]  = cspice_bodvrd('MOON', 'RADII', 3)*1E3;
% % R_M  = radii(1);

%% MOON VALUES FROM GRMGM 12000 FILE
R_M      = 1.7380000000000000e+06;   % [m]
GM_moon  = 4.9028001224453001e+12;   % [m^3/ s^2]
GM_sigma = 6.4536052689015518e-24;   % [m^3/ s^2]

%% RMS VALUE COEFFICIENT
[RMS_SH_coeff] = computeRMS_coeffErr(n_max, Nc, Ns, SH_coeff, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_error] = computeRMS_coeffErr(n_max, Nc, Ns, SH_clone, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_uncrt] = computeRMS_coeffErr(n_max, Nc, Ns, SH_uncrt, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

figure()
semilogy(2:n_max, RMS_SH_coeff(2:end), 'LineWidth', 2, 'Color', 'b'); 
grid on; hold on;
semilogy(2:n_max, 3.*RMS_SH_uncrt(2:end), 'LineWidth', 2, 'Color', 'b', ...
    'LineStyle', '--'); 
semilogy(2:n_max, RMS_SH_error(2:end), 'LineWidth', 2, 'Color', 'r', ...
    'LineStyle', '-'); 
xlabel('Degree'); title('RMS value GRMG1200');
legend('RMS coefficients', 'RMS \sigma'); 

%% LRO SPICE trajectory
utc_start = '2012-03-01 00:00:00';
utc_stop  = '2012-03-04 00:00:00';

et0       = cspice_str2et(utc_start);
et1       = cspice_str2et(utc_stop);
time_sec  = et1 - et0;

frec      = 1 /(10*60); % [Hz]
N         = round(time_sec * frec);   % number of samples
N         = 10;

et        = linspace(et0, et1, N);

tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sc_SPICE(1:3, :) = sc_SPICE(1:3, :).*1E3;         % [m]
sc_SPICE(4:6, :) = sc_SPICE(4:6, :).*1E3;         % [m/s]

% convert time to UTC
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

figure()
plot3(sc_SPICE(1, :), sc_SPICE(2, :), sc_SPICE(3, :), 'LineWidth', 2)
hold on; grid on; title('LRO SPICE orbit'); axis equal;

figure()
plot(tUTC, (vecnorm(sc_SPICE(1:3, :)) - R_M)./1E3, 'LineWidth', 2);
title('Orbit altitude');ylabel('[km]')

disp('Computing rotation matrices')
for j = 1:N
    frame_to   = 'J2000';
    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et(j));

    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000_mat(minInd:maxInd, :) = J2000_MOON';
end
%% GG measurements (true + perturbed)
Monte_Carlo = 1; % number of monte carlo realizations
n_max       = 1199; [Nc, Ns, Ncs]    = count_num_coeff(1200);

Y_true      = nan(9, N, Monte_Carlo);
Y_nom       = nan(9, N, Monte_Carlo);

% trajectory points
x = sc_SPICE(1, :); y = sc_SPICE(2, :); z = sc_SPICE(3, :);

scale_dist = 1E3;
x = x / scale_dist; y = y /scale_dist; z = z /scale_dist;
R_M = R_M / scale_dist;
GM_moon = GM_moon / (scale_dist^3);
GM_sigma = GM_sigma/ (scale_dist^3);

disp('Computing reference gravity field');
[Cnm, Snm] = list2mat(1200, Nc, Ns, SH_coeff);
for j = 1:N
    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);

    r = [x(j);y(j);z(j)]; v = zeros(3, 1);                            
    [Y_J2000, ~] = gradiometer_meas(et(j) ,...
        [GM_moon, R_M, n_max, normalized],...
        BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);

    % rotate to ENU coordinates
    ENU_ECEF = ecef2enu(BODYMOON_J2000 * r);
    T_J2000  = reshape(Y_J2000, [3,3]);
    T_ECEF   = BODYMOON_J2000 * T_J2000 * BODYMOON_J2000';
    T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';

    Y_ENU    = reshape(T_ENU', [9, 1]);

    Y_nom(:, j, :) = Y_ENU.*ones(9, Monte_Carlo);
end

disp('Computing MC simulation')
for k = 1:Monte_Carlo
    disp("computing ... "  + string(k/Monte_Carlo * 100) + "%")
    err = normrnd(0, SH_uncrt(2:end));
    SH  = SH_coeff + [0;SH_clone(2:end)];
    GM_moon_2 = GM_moon + normrnd(0, GM_sigma);

    [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH);
    tic;
    for j = 1:N
        maxInd = 3 *j; minInd = maxInd - 2;
        BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);
    
        r = [x(j);y(j);z(j)];
        v = zeros(3, 1);                            
        [Y_J2000, ~] = gradiometer_meas(et(j) ,...
            [GM_moon_2, R_M, n_max, normalized],...
            BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);

        % rotate to ENU coordinates
        ENU_ECEF = ecef2enu(BODYMOON_J2000 * r);
        T_J2000  = reshape(Y_J2000, [3,3]);
        T_ECEF   = BODYMOON_J2000 * T_J2000 * BODYMOON_J2000';
        T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
        Y_ENU    = reshape(T_ENU', [9, 1]);
    

        Y_true(:, j, k) = Y_ENU;
    end
    sim_t = toc;
    if(k == 1)
        disp('  Estimated time... ' + ...
        string(sim_t * Monte_Carlo / 60) + ' min')
    end
end

disturbance = (Y_true - Y_nom)./1E-12; % [milli-Eotvos]

% stats
RMS_distubance_xx = rms(squeeze(disturbance(1, :, :)), 2);
RMS_distubance_xy = rms(squeeze(disturbance(2, :, :)), 2);
RMS_distubance_xz = rms(squeeze(disturbance(3, :, :)), 2);
RMS_distubance_yy = rms(squeeze(disturbance(5, :, :)), 2);
RMS_distubance_yz = rms(squeeze(disturbance(6, :, :)), 2);
RMS_distubance_zz = rms(squeeze(disturbance(9, :, :)), 2);

RMS_distubance    = [RMS_distubance_xx, RMS_distubance_xy,RMS_distubance_xz ...
    RMS_distubance_yy,RMS_distubance_yz, RMS_distubance_zz];
disp('Time mean for RMS disturbance (milli-Eotvos):');
disp(mean(RMS_distubance));

% plot resuts
figure()
tt = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", 'zz'];
for k = 1:9
    subplot(3, 3, k)
    plot(tUTC, squeeze(disturbance(k, :, :)), '.', 'Color', 'g');
    title(tt(k)); grid on; ylabel("[mE]")
end
