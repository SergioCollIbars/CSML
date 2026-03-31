clear; clc; close all;
format long g;
addpath('../../data/'); addpath(genpath('../../functions'));
set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')
%%          ATITTUDE ERROR BUDGET LRO
% Description: Compute the impact of gravity field coefficient uncertainty
% in the LRO orbit. (error budget)
% Author: Sergio Coll-Ibars
% Date: 01/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRGM1200 gravity field 
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
R_M = file(2)*1E3; normalized = file(3);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);

%% LRO SPICE trajectory
utc_start = '2012-03-20 00:00:00';
utc_stop  = '2012-03-21 00:00:00';
N         = 1000;                         % number of samples
[GM] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Moon [km^3/s^2]
GM_moon = GM * 1E9;                       % [m^3/s^2]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

dt = et(2) - et(1);

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

% rotation matrices
BODYMOON_J2000_mat = nan(3 * N , 3);
BODYSC_J2000_mat   = nan(3 * N, 3);
for j = 1:N
    frame_to   = 'J2000';
    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et(j));

    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000_mat(minInd:maxInd, :) = J2000_MOON';

    J2000_RTN   = RTN2ECI(sc_SPICE(1:3, j), sc_SPICE(4:6, j));
    BODYSC_J2000_mat(minInd:maxInd, :) = J2000_RTN';
end

% compute angular velocity & error
frec = 1;
[angVel_vec] = angularVel_from_DCM(BODYSC_J2000_mat, dt);
radians_per_sec = 0.006 * (pi/180) / 3600 / sqrt(1/frec);

% plot angular velocity
figure();
semilogy(tUTC(2:end-1), abs(angVel_vec(:, 2:end-1)), 'LineWidth', 2);
grid on; ylabel('[rad/s]');
legend('\omega_R', '\omega_T', '\omega_N');
title('Angular velocity in Body frame coordinates');

% plot orbit altitude
figure();
plot(tUTC, (vecnorm(sc_SPICE(1:3, :)) - R_M)./1E3, 'LineWidth', 2);
grid on; ylabel('[km]');
title('Orbit altitude, R = ' + string(R_M./1E3) + ' Km');

%% GG measurements (true + perturbed)
Monte_Carlo = 100; % number of monte carlo realizations
n_max       = 100; [Nc, Ns, Ncs]    = count_num_coeff(n_max);

[Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);
disturbance = nan(6, N, Monte_Carlo);
for k = 1:Monte_Carlo
    err_angVel = normrnd(0, radians_per_sec, [3, N]);

    for j = 2:N-1
        [W1] = dyad_operator(angVel_vec(:, j) + err_angVel(:, j));
        [W2] = dyad_operator(angVel_vec(:, j));
        
        dW = W1*W1 - W2*W2;
        disturbance(:, j, k) = [dW(1,1);dW(1,2);dW(1,3);dW(2,2);dW(2,3);...
            dW(3,3)]./1E-12;    % [milli-Eotvos]
    end
end

RMS_distubance_xx = rms(squeeze(disturbance(1, :, :)), 2, 'omitnan');
RMS_distubance_xy = rms(squeeze(disturbance(2, :, :)), 2, 'omitnan');
RMS_distubance_xz = rms(squeeze(disturbance(3, :, :)), 2, 'omitnan');
RMS_distubance_yy = rms(squeeze(disturbance(4, :, :)), 2, 'omitnan');
RMS_distubance_yz = rms(squeeze(disturbance(5, :, :)), 2, 'omitnan');
RMS_distubance_zz = rms(squeeze(disturbance(6, :, :)), 2, 'omitnan');

RMS_distubance    = [RMS_distubance_xx, RMS_distubance_xy,RMS_distubance_xz ...
    RMS_distubance_yy,RMS_distubance_yz, RMS_distubance_zz];
disp('Time mean for RMS disturbance (milli-Eotvos):');
disp(mean(RMS_distubance, 'omitnan'));

%% FUNCTIONS
function [M] = dyad_operator(x)
    % construct dyad operator
    M = [0, -x(3), x(2);...
        x(3), 0, -x(1);...
        -x(2), x(1), 0];
end