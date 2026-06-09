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
n_max      = 100; 

% S/C mass and cross-area
A = 2.8;      % [m^2]
M = 220;      % [kg]

%% LRO SPICE trajectory
utc_start = '2015-03-20 00:00:00 UTC';
utc_stop  = '2015-03-20 12:00:00 UTC';
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
% % tgt       = 'GRAIL-B';
tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

r_sc = sc_SPICE(1:3, :).*1E3;         % [m]
v_sc = sc_SPICE(4:6, :).*1E3;         % [m/s]

% plot altitude
figure();
H = (vecnorm(r_sc) - 1738E3)./1E3;
plot(tUTC, H, 'LineWidth', 2);
grid on; ylabel('[Km]');


% Compute different forces norm
FORCES   = nan(8, length(et));
GRADIENT = nan(4, length(et));

for k = 1:length(et)

    t = et(k);

    % ================================================================
    % State and geometry
    % ================================================================

    % Spacecraft state relative to Moon, expressed in J2000
    r_sc_M = r_sc(:, k);       % Moon -> spacecraft [m]
    v_sc_M = v_sc(:, k);       % spacecraft wrt Moon [m/s]

    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';

    % Moon -> Earth, J2000
    [Estate, ~] = cspice_spkezr('EARTH', t, ref, abcorr, observer);
    r_E_M = Estate(1:3) * 1e3;       % [m]

    % Moon -> Sun, J2000
    [Sstate, ~] = cspice_spkezr('SUN', t, ref, abcorr, observer);
    r_S_M = Sstate(1:3) * 1e3;       % [m]

    % Relative vectors
    r_SC_E = r_sc_M - r_E_M;         % Earth -> spacecraft [m]
    r_M_E  = -r_E_M;                 % Earth -> Moon [m]

    r_SC_S = r_sc_M - r_S_M;         % Sun -> spacecraft [m]
    r_M_S  = -r_S_M;                 % Sun -> Moon [m]

    % Moon-centered spacecraft vector
    r_SC_M = r_sc_M;                 % Moon -> spacecraft [m]

    % ================================================================
    % Frame transformations
    % ================================================================

    J2000_EARTH = cspice_pxform('IAU_EARTH', 'J2000', t);
    J2000_MOON  = cspice_pxform('MOON_PA',   'J2000', t);

    % ================================================================
    % Shadow and SRP
    % ================================================================

    % shadow_function expects:
    %   r_sun = Moon -> Sun
    %   r_sc  = Moon -> spacecraft
    F = shadow_function(R_Sun(1), R_M, r_S_M, r_SC_M);

    % SRP function assumed to expect Sun -> spacecraft vector
    [aSRP, daSRP_dr, ~] = SRP(r_SC_S, 1.75, M, A);

    % ================================================================
    % Earth gravity: two-body and oblateness
    % ================================================================

    Cmat_E = Cnm_list{1};
    Smat_E = Snm_list{1};

    % Earth acceleration on spacecraft, degree 2
    [~, a_E_SC_EF, ddU1_EF] = potentialGradient_nm( ...
        Cmat_E, Smat_E, 2, ...
        J2000_EARTH' * r_SC_E, R_E, GM_E, normalized);

    % Earth central acceleration on spacecraft
    [~, a_E_SC_0_EF, ~] = potentialGradient_nm( ...
        Cmat_E, Smat_E, 0, ...
        J2000_EARTH' * r_SC_E, R_E, GM_E, normalized);

    % Earth central acceleration on Moon
    [~, a_E_M_0_EF, ~] = potentialGradient_nm( ...
        Cmat_E, Smat_E, 0, ...
        J2000_EARTH' * r_M_E, R_E, GM_E, normalized);

    % Rotate all Earth accelerations back to J2000
    a_E_SC      = J2000_EARTH * a_E_SC_EF;
    a_E_SC_0    = J2000_EARTH * a_E_SC_0_EF;
    a_E_M_0     = J2000_EARTH * a_E_M_0_EF;
    ddU1        = J2000_EARTH * ddU1_EF * J2000_EARTH';

    % Earth differential two-body acceleration in Moon-centered dynamics
    a_Earth_2body = a_E_SC_0 - a_E_M_0;

    % Earth non-central/oblateness contribution on spacecraft
    a_Earth_oblate = a_E_SC - a_E_SC_0;

    % ================================================================
    % Lunar gravity: central and non-central
    % ================================================================

    Cmat_M = Cnm_list{2};
    Smat_M = Snm_list{2};

    % Full Moon gravity to n_max
    [~, a_M_SC_MF, ddU2_MF] = potentialGradient_nm( ...
        Cmat_M, Smat_M, n_max, ...
        J2000_MOON' * r_SC_M, R_M, GM_M, normalized);

    % Moon central gravity
    [~, a_M_SC_0_MF, ~] = potentialGradient_nm( ...
        Cmat_M, Smat_M, 0, ...
        J2000_MOON' * r_SC_M, R_M, GM_M, normalized);

    % Non-central lunar gravity
    a_Moon_oblate_MF = a_M_SC_MF - a_M_SC_0_MF;

    % Rotate back to J2000
    a_Moon_central = J2000_MOON * a_M_SC_0_MF;
    a_Moon_oblate  = J2000_MOON * a_Moon_oblate_MF;
    ddU2           = J2000_MOON * ddU2_MF * J2000_MOON';

    % ================================================================
    % Sun third-body gravity
    % ================================================================

    % Direct point-mass implementation.
    % This avoids sign ambiguity from potentialGradient_nm.
    a_Sun_SC = -GM_S * r_SC_S / norm(r_SC_S)^3;   % Sun acceleration on SC
    a_Sun_M  = -GM_S * r_M_S  / norm(r_M_S)^3;    % Sun acceleration on Moon

    % Differential Sun acceleration in Moon-centered dynamics
    a_Sun_3body = a_Sun_SC - a_Sun_M;

    % Approximate Sun gravity-gradient tensor at spacecraft, J2000
    I3 = eye(3);
    rhoS = norm(r_SC_S);
    ddU3 = GM_S * (3 * (r_SC_S * r_SC_S.') / rhoS^5 - I3 / rhoS^3);

    % ================================================================
    % Relativity
    % ================================================================

    % Simple Moon Schwarzschild correction
    a_rel = relativistic_accel(r_SC_M, v_sc_M, GM_M);

    % ================================================================
    % Lunar albedo
    % ================================================================

    % lunar_albedo_accel expects:
    %   r_sc_M = Moon -> spacecraft
    %   r_S_M  = Moon -> Sun
    a_alb = lunar_albedo_accel(r_SC_M, r_S_M, 1.3, A, M, 0.12);

    % ================================================================
    % Store acceleration magnitudes
    % ================================================================

    FORCES(1, k) = norm(a_Earth_2body);
    FORCES(2, k) = norm(a_Moon_central);
    FORCES(3, k) = norm(a_Sun_3body);
    FORCES(4, k) = norm(a_Moon_oblate);
    FORCES(5, k) = norm(a_Earth_oblate);
    FORCES(6, k) = norm(a_rel);
    FORCES(7, k) = norm(F * a_alb);
    FORCES(8, k) = norm(F * aSRP);

    % ================================================================
    % Store gradient magnitudes in Eotvos
    % ================================================================

    GRADIENT(1, k) = norm(ddU1, 'fro') / 1e-12;
    GRADIENT(2, k) = norm(ddU2, 'fro') / 1e-12;
    GRADIENT(3, k) = norm(ddU3, 'fro') / 1e-12;
    GRADIENT(4, k) = norm(F * daSRP_dr, 'fro') / 1e-12;

end
% plots
figure()
semilogy(tUTC, abs(FORCES(1:7, :))./1E3, 'LineWidth', 2);hold on;
semilogy(tUTC, abs(FORCES(8, :))./1E3,'LineStyle', 'none', 'Marker','.');
ylabel('km/s^2'); grid on; 
title('LRO accelerations');
legend('Earth Two-body', 'Central-Moon', 'Sun gravity', 'Lunar Oblate', ...
    'Earth Oblate', 'Relativity', 'Lunar albedo','SRP');


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

function a_rel = relativistic_accel(r, v, GM)

    % Relativistic Schwarzschild correction
    %
    % Inputs:
    %   r  : spacecraft position relative to central body [m]
    %   v  : spacecraft velocity relative to central body [m/s]
    %   GM : central body gravitational parameter [m^3/s^2]
    %
    % Output:
    %   a_rel : relativistic acceleration [m/s^2]

    c = 299792458;   % speed of light [m/s]

    rnorm = norm(r);
    v2    = dot(v, v);
    rv    = dot(r, v);

    a_rel = GM/(c^2 * rnorm^3) * ...
        ( (4*GM/rnorm - v2)*r + 4*rv*v );
end

function a_alb = lunar_albedo_accel(r_sc_M, r_S_M, CR, Area, mass, albedo)
    % Lunar albedo acceleration, simple cannonball approximation
    %
    % Inputs:
    %   r_sc_M  : Moon -> spacecraft vector [m], J2000
    %   r_S_M   : Moon -> Sun vector [m], J2000
    %   CR      : reflectivity coefficient [-]
    %   Area    : effective cross-sectional area [m^2]
    %   mass    : spacecraft mass [kg]
    %   albedo  : lunar albedo coefficient [-]
    %
    % Output:
    %   a_alb   : lunar albedo acceleration [m/s^2], J2000
    
    c  = 299792458;        % [m/s]
    S0 = 1361;             % solar constant at 1 AU [W/m^2]
    AU = 149597870700;     % [m]
    R_M = 1738e3;          % mean lunar radius [m]
    
    r = norm(r_sc_M);
    
    % Correct solar flux for Moon-Sun distance
    r_MS = norm(r_S_M);
    S = S0 * (AU/r_MS)^2;
    
    % Unit vectors
    u_M_to_SC = r_sc_M / norm(r_sc_M);
    u_M_to_S  = r_S_M  / norm(r_S_M);
    
    % Phase-angle factor: illuminated fraction visible from spacecraft
    cos_alpha = dot(u_M_to_SC, u_M_to_S);
    
    % Keep only visible illuminated contribution
    phase = max(0, (1 + cos_alpha)/2);
    
    % Radiation pressure from reflected lunar light
    P_alb = (S/c) * albedo * (R_M/r)^2 * phase;
    
    % Acceleration direction: away from Moon
    a_alb = CR * Area/mass * P_alb * u_M_to_SC;
end