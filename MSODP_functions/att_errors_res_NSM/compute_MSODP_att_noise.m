clear; clc;
close all;
format long g;

%% COMPUTE ATTITUDE NOISE from the nominal vs true observation attitude 
% Files generated with MSODP
% Author: Sergio Coll-Ibars
% Date: 03/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILE_noise_free = "/Users/sergiocollibars/Documents/att_files/" + ...
    "grace.2008-08-01_A_61_nonoise.att.pre";
att_noise_free  = read_att_to_rpy(FILE_noise_free);
DCM_noise_free  = read_att_to_DCM(FILE_noise_free);

FILE_noise = "/Users/sergiocollibars/Documents/att_files/" + ...
    "grace.2008-08-01_A_61_mode3.att.pre";
att_noise = read_att_to_rpy(FILE_noise);
DCM_noise = read_att_to_DCM(FILE_noise);

Nt      = length(DCM_noise(1, 1, :));
DCM_err = nan(3,3, Nt);
for k = 1:Nt
    DCM_err(:,:, k) =  DCM_noise_free(:,:,k)' * DCM_noise(:, :, k); % convetion in MSODP: R_actual = R_nom * R_noise
end

err_rpy = nan(3, Nt); % yaw-pitch-roll (3-2-1)
for k =1:Nt
    C = DCM_err(:,:, k);
    y = atan2(C(1,2), C(1,1)); % yaw
    p = -asin(C(3,1));         % pitch
    r = atan2(C(2,3), C(3,3)); % roll

    err_rpy(:, k) = [r;p;y];
end

figure;
plot(att_noise.t, err_rpy(1, :))
xlabel('Time since first sample [s]')
ylabel('Roll [rad]')
grid on

figure;
plot(att_noise.t, err_rpy(2, :))
xlabel('Time since first sample [s]')
ylabel('Pitch [rad]')
grid on

figure;
plot(att_noise.t, err_rpy(3, :))
xlabel('Time since first sample [s]')
ylabel('Yaw [rad]')
grid on

figure;
plot(att_noise.t,  att_noise_free.roll)
xlabel('Time since first sample [s]')
ylabel('Roll [rad]')
grid on

figure;
plot(att_noise.t,att_noise_free.pitch)
xlabel('Time since first sample [s]')
ylabel('Pitch [rad]')
grid on

figure;
plot(att_noise.t, att_noise_free.yaw)
xlabel('Time since first sample [s]')
ylabel('Yaw [rad]')
grid on

RMS_roll  = rms(err_rpy(1, :));
RMS_pitch = rms(err_rpy(2, :));
RMS_yaw   = rms(err_rpy(3, :));