clear; clc; close all;
format long g;

addpath('../data/'); addpath(genpath("../functions/"));
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
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
n_max            = file(1); R_M = file(2); normalized = file(3);
[Nc, Ns, ~]    = count_num_coeff(n_max);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);

% compute RMS value
[RMS_SH_coeff] = computeRMS_coeffErr(n_max, Nc, Ns, SH_coeff, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_uncrt] = computeRMS_coeffErr(n_max, Nc, Ns, SH_uncrt, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

figure()
semilogy(2:n_max, RMS_SH_coeff(2:end), 'LineWidth', 2, 'Color', 'b'); 
grid on; hold on;
semilogy(2:n_max, RMS_SH_uncrt(2:end), 'LineWidth', 2, 'Color', 'b', ...
    'LineStyle', '--'); 
xlabel('Degree'); title('RMS value GRMG1200');
legend('RMS coefficients', 'RMS \sigma'); 

%% LRO SPICE trajectory
utc_start = '2025-03-16 04:30:00';
utc_stop  = '2025-03-17 00:00:00';
N         = 200;                         % number of samples
[GM] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Moon [km^3/s^2]
GM_moon = GM * 1E9;                       % [m^3/s^2]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

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

for j = 1:N
    frame_to   = 'J2000';
    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et(j));

    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000_mat(minInd:maxInd, :) = J2000_MOON';
end
%% GG measurements (true + perturbed)
Monte_Carlo = 100; % number of monte carlo realizations
n_max       = 100; [Nc, Ns, Ncs]    = count_num_coeff(n_max);

Y_true      = nan(9, N, Monte_Carlo);
Y_nom       = nan(9, N, Monte_Carlo);

[Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);
for j = 1:N
    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);

    r = sc_SPICE(1:3, j); v = sc_SPICE(4:6, j);                            
    [Y_J2000, ~] = gradiometer_meas(et(j) ,...
        [GM_moon, R_M, n_max, normalized],...
        BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);
    Y_nom(:, j, :) = Y_J2000.*ones(9, Monte_Carlo);
end

for k = 1:Monte_Carlo
    err = normrnd(0, SH_uncrt(2:end));
    SH  = SH_coeff + [0;err];

    [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH);
    for j = 1:N
        maxInd = 3 *j; minInd = maxInd - 2;
        BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);
    
        r = sc_SPICE(1:3, j); v = sc_SPICE(4:6, j);                            
        [Y_J2000, ~] = gradiometer_meas(et(j) ,...
            [GM_moon, R_M, n_max, normalized],...
            BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);
        Y_true(:, j, k) = Y_J2000;
    end
end

disturbance = (Y_true - Y_nom)./1E-12; % [milli-Eotvos]

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