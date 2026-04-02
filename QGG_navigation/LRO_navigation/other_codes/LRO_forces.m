clear;
clc;
close all;

addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data"))
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/functions"))

set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');
%% FORCE QUANTIFICATION IN LRO  
% Plot the forces from different sources and the tensor for the LRO;
% Date: 04/02/2026
% Author: Sergio Coll-Ibars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start Simulation
metaData_path = "Metadata.txt";
mtd = readParams("data/"+metaData_path);
fld = 1;
[planetParams, Cnm_list, Snm_list, ...
    state0, t_range]              = loadUniverse(mtd.folder{fld});


GM_E = planetParams(1); GM_M = planetParams(2);
R_E  = planetParams(3); R_M  = planetParams(4);
[GM_S] = cspice_bodvrd('SUN', 'GM', 1)*1E9; 

% gravity computation params
normalized = planetParams(6);
n_max      = 0; 

% S/C mass and cross-area
A = 2.8;      % [m^2]
M = 220;      % [kg]

%% LRO SPICE trajectory
utc_start = '2012-03-04 00:00:00';
utc_stop  = '2012-03-06 00:00:00';
N         = 20000;                           % number of samples
[R_Sun] = cspice_bodvrd('SUN', 'RADII', 3).*1E3;      % Get R for the Sun  [m]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

dt = et(2) - et(1);

% convert time to UTC
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

% S/C trajectroy
tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

r_sc = sc_SPICE(1:3, :).*1E3;         % [m]
v_sc = sc_SPICE(4:6, :).*1E3;         % [m/s]


% Compute different forces norm
FORCES = nan(4, length(et)); GRADIENT = nan(4, length(et));
for k = 1:length(et)
    x = r_sc(:, k);

    % compute Earth position. ref: J2000
    target = 'EARTH';
    t = et(k);                         % Convert UTC time to ephemeris time

    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';                 % Set the observer to the MOON

    [Estate, ~] = cspice_spkezr(target, t, ref, abcorr, observer);     % [Km & Km/s]
    Epos  = Estate(1:3)*1E3;                                            % [m]
    r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

    % compute Moon position. ref: MOON
    target = 'MOON';
    [Mstate, ~] = cspice_spkezr(target, t, ref, abcorr, observer);     % [Km & Km/s]
    Mpos  = Mstate(1:3)*1E3;                                            % [m]
    r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

    % compute Sun position. ref: MOON
    target = 'SUN';
    [Sstate, ~] = cspice_spkezr(target, t, ref, abcorr, observer);     % [Km & Km/s]
    Spos = Sstate(1:3)*1E3;                                             % [m]
    r3 = [x(1)-Spos(1);x(2)-Spos(2);x(3)-Spos(3)];                      % SC-Sun
    
    % relative vector from Earth & Sun to Moon
    r_ME = Mpos - Epos;
    r_MS = Mpos - Spos;

    % compute orientation
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, t);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, t);

    % SRP acceleration
    F = shadow_function(R_Sun(1), R_M, Spos, r2);
    [aSRP, daSRP_dr, ~] = SRP(r3, 2, M, A);

    % compute gravity acceleration
    Cmat_E = Cnm_list{1};
    Smat_E = Snm_list{1};
    [~, dU1, ddU1] = potentialGradient_nm(Cmat_E, Smat_E, 0, ...
                                                J2000_EARTH'*r1, R_E, GM_E, ...
                                                normalized);

    [~, dU1_T, ~] = potentialGradient_nm(Cmat_E, Smat_E, 0, ...
                                                J2000_EARTH'*r_ME, R_E, GM_E, ...
                                                normalized);

    Cmat_M = Cnm_list{2};
    Smat_M = Snm_list{2};

    [~, dU2, ddU2] = potentialGradient_nm(Cmat_M, Smat_M, n_max, ...
                                                J2000_MOON'*r2, R_M, GM_M, ...
                                                normalized);


    [~, dU3, ddU3] = potentialGradient_nm(Cmat_M, Smat_M, 0, ...
                                                r3, R_Sun(1), GM_S, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    dU1  = J2000_EARTH  * dU1;
    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';

    dU2  = J2000_MOON  * dU2;
    ddU2 = J2000_MOON  * ddU2  * J2000_MOON'; 

    % Tidial acceleration
    a_tidial_E = - J2000_EARTH  * dU1_T;
    a_tidial_S = GM_S * r_MS / (vecnorm(r_MS)^3);
    
    % total acceleration
    FORCES(1, k) = vecnorm(dU1 + a_tidial_E);
    FORCES(2, k) = vecnorm(dU2);
    FORCES(3, k) = vecnorm(dU3 + a_tidial_S);
    FORCES(4, k) = vecnorm(F * aSRP);


    GRADIENT(1, k) = norm(ddU1, 'fro')./1E-12;
    GRADIENT(2, k) = norm(ddU2, 'fro')./1E-12;
    GRADIENT(3, k) = norm(ddU3, 'fro')./1E-12;
    GRADIENT(4, k) = norm(F * daSRP_dr, 'fro')./1E-12;
end

% plots
figure()
semilogy(tUTC, abs(FORCES)./1E3, 'LineWidth', 2);
ylabel('km/s^2'); grid on;
title('LRO accelerations');
legend('Earth gravity', 'Lunar gravity', 'Sun gravity', 'SRP');


figure()
semilogy(tUTC, GRADIENT, 'LineWidth', 2);
ylabel('milli-Eotvos'); grid on;
title('LRO gradient');
legend('Earth gradient', 'Lunar gradient', 'Sun gradient', ...
    'SRP gradient');



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