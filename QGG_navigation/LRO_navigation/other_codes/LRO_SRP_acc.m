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

% S/C Mass and cross-Area
A = 1.09*0.95;      % [m^2]
M = 202.4;          % [kg]

%% LRO SPICE trajectory
utc_start = '2012-03-19 00:00:00';
utc_stop  = '2012-03-19 04:00:00';
N         = 2000;                         % number of samples
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
grid on; title('Shadow Function'); grid on;

% compute SRP and SRP tensor acceleration
F_SRP = nan(3, length(tUTC)); T_SRP = nan(9, length(tUTC));
for k = 1:length(tUTC)
    F = shadow_function(R_Sun(1),R_Moon(1),sun_SPICE(1:3, k),...
        sc_SPICE(1:3, k));

     % relative vector from Sun to S/C
    r3 = sc_SPICE(1:3, k) - sun_SPICE(1:3, k);
    [aSRP, daSRP_dr, ~] = SRP(r3, 1, M, A);

    F_SRP(:, k) = F * aSRP;
    T_SRP(:, k) = reshape(F*daSRP_dr, [9, 1]);
end

figure()
plot(tUTC, F_SRP, 'LineWidth', 2); grid on;
title('SRP acceleration'); ylabel('m/s^2');

figure()
plot(tUTC, norm(T_SRP)./1E-12, 'LineWidth', 2); grid on;
title('SRP gradient'); ylabel('milli-Eotvos');

%% AUXILIARY FUNCTIONS
function [aSRP, daSRP_dr, daSRP_dEta] = SRP(rs, eta, m, A)
    %%                   SRP ACCELERATION FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 01/20/2023
    %
    %   Description: This function computes the SRP acceleration at certain
    %   time.
    %
    %   Input:
    %       rs: SC position to sun inerital frame. ECI frame
    %       eta: scale factor
    %       A: S/C reference area [-]
    %       m: S/C reference mass [kg]
    %
    %   Output: 
    %       aSRP: SRP acceleration. ACI frame
    %       daSRP_dr: partial respect to inertial position
    %       daSRP_dEta: partial respect to eta
    % --------------------------------------------------------------------%

    % SRP parameters
    Gamma = 5.67E-8;                % [Kg/K^4]
    c     = 2.99792458E8;               % [m/s^2] 
    Rs    = 6.96E8 ;                   % [m]
    Ts    = 5778;                      % [K]
    Cr    = 1;

    % SRP acc
    rsn = vecnorm(rs);
    P = (Gamma * Rs^2 * Ts^4 )/ (c * m);

    aSRP = eta * ( P * Cr * A) * rs / rsn^3;
    
    % SPR partial respect to inertial position
    daSRP_dr = eta * Gamma * Ts^4 * Rs^2 * Cr * A / (m * c * rsn^3) * ...
        (eye(3) - 3 * (rs * rs') / (rsn^2));

    % SRP partial respect to eta
    daSRP_dEta = ( P * Cr * A) * rs / rsn^3;
end