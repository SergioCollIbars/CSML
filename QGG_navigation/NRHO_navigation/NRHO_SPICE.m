clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")
%%                 QGG SPICE NRHO ORBIT
% Description: test the ODE 113 integration accuracy vs the SPICE generated
% orbit. Tested orbit: Lunar gateway NRHO. 

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

ET0  = cspice_str2et({'2022-11-13 00:00:00 TDB'}); 
ET1  = cspice_str2et({'2023-05-18 00:00:00 TDB'});

frec = 1 / 10;  % [Hz]
Nt   = round((ET1 - ET0) * frec);
time = linspace(ET0, ET1, Nt);

[sc_SPICE, ~] = cspice_spkezr('-60000', time, 'J2000', 'NONE', '3');

Xsc_SPICE = sc_SPICE(1, :);  Ysc_SPICE = sc_SPICE(2, :);
Zsc_SPICE = sc_SPICE(3, :);

VXsc_SPICE = sc_SPICE(4, :);  VYsc_SPICE = sc_SPICE(5, :);
VZsc_SPICE = sc_SPICE(6, :);

% convert time to date
jd = 2451545 + time / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', 'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';
date = humanReadableTime;

% plot SPICE trajectory
figure()
plot3(sc_SPICE(1, :) , sc_SPICE(2, :) , sc_SPICE(3, :), 'LineWidth', 2)
axis equal;
grid on; title('SPICE 3D trajectory J2000 frame');

figure()
plot(date, vecnorm(sc_SPICE(1:3, :)), 'LineWidth', 2);
title('S/C orbit radius norm');

% save trajectory
traj_pos = [Xsc_SPICE;Ysc_SPICE;Zsc_SPICE];                                % [km]
traj_vel = [VXsc_SPICE;VYsc_SPICE;VZsc_SPICE];                             % [km/s]
t        = time;                                                           % Epoch [sec]
traj = [traj_pos;traj_vel;t];
name = "CAPSTONE_11_13_2022_to_05_18_2023.mat";
save(name, "traj")

% close SPICE
cspice_kclear