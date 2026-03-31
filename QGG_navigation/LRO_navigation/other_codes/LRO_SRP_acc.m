clear; clc; close all;
format long g;
addpath('/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data'); 
addpath(genpath('/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/functions'));
set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')
%%          SRP ACCELERATION LRO
% Description: Compute the SRP acceleration accounting for shadows
% Author: Sergio Coll-Ibars
% Date: 01/27/2026
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
utc_start = '2012-03-19 04:30:00';
utc_stop  = '2012-03-19 06:00:00';
N         = 20000;                         % number of samples
[GM] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Moon [km^3/s^2]
GM_moon = GM * 1E9;                       % [m^3/s^2]

[R_Moon] = cspice_bodvrd('MOON', 'RADII', 3).*1E3;    % Get R for the Moon [m]
[R_Sun] = cspice_bodvrd('SUN', 'RADII', 3).*1E3;      % Get R for the Sun  [m]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

dt = et(2) - et(1);

% S/C trajectroy
tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sc_SPICE(1:3, :) = sc_SPICE(1:3, :).*1E3;         % [m]
sc_SPICE(4:6, :) = sc_SPICE(4:6, :).*1E3;         % [m/s]

% Sun trajectory
tgt       = 'SUN';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sun_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sun_SPICE(1:3, :) = sun_SPICE(1:3, :).*1E3;         % [m]
sun_SPICE(4:6, :) = sun_SPICE(4:6, :).*1E3;         % [m/s]

% convert time to UTC
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");


% compute shadow function
F = nan(1, length(tUTC));
for k = 1:length(tUTC)
    F(k) = shadow_function(R_Sun(1),R_Moon(1),sun_SPICE(1:3, k),...
        sc_SPICE(1:3, k));
end

figure()
plot(tUTC, F, 'LineWidth', 2);
grid on; title('Shadow Function');