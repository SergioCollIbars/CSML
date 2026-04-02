clear;
clc;
close all;

addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/data"))
addpath(genpath("/Users/sergiocollibars/Desktop/CSML/QGG_navigation/LRO_navigation/functions"))


set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');
%% SRP FORCE IN EOM INTEGRATOR
% Test the inclusiom of a SRP model to the Numerical integrator for Low
% Lunar Orbits.
% Date: 04/01/2026
% Author: Sergio Coll Ibars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load meta data
metaData_path = "Metadata.txt";
mtd = readParams("data/"+metaData_path);
fld = 1;

% Spacecraft Mass & Area
A = 1.09*0.95;      % [m^2]
M = 202.4;          % [kg]

 % start Simulation
[planetParams, Cnm_list, Snm_list, ...
    state0, t_range]              = loadUniverse(mtd.folder{fld});
[instrumentParams_GG]             = loadInstrument_GG(mtd.folder{fld});
[instrumentParams_ST]             = loadInstrument_ST(mtd.folder{fld});

% integrate Trajectory
disp('Simulating trajectory ...')
tmin = t_range(1);
tmax = t_range(2);
t = linspace(tmin, tmax, (tmax-tmin) * instrumentParams_GG(1, 5));

options = odeset('RelTol',1e-13,'AbsTol',1e-13); Nx = 6;
PHI0    = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
X0      = state0;
[time, state] = ode113(@(t, x) EOM_LRO_EPHEM_SRP_TEST(t, x, planetParams, ...
    Cnm_list, Snm_list, A, M), t, [X0;PHI0], ...
    options);

disp('  DONE ...')

% Plot results
figure(); ll = ["X", "Y", "Z"];
for j = 1:3
    subplot(1, 3, j);
    plot(time - time(1), state(:, j)./1E3, 'LineWidth', 2); grid on;
    ylabel('[km]'); title('J200 Position ' + ll(j));
    xlabel('Epoch [s]');
end

figure();
for j = 1:3
    subplot(1, 3, j);
    plot(time - time(1), state(:, j+3)./1E3, 'LineWidth', 2); grid on;
    ylabel('[km/s]'); title('J200 Velocity ' + ll(j));
    xlabel('Epoch [s]');
end
figure()
plot(time - time(1), (vecnorm(state(:, 1:3)') - 1738E3)./1E3, 'LineWidth', 2); 
grid on; xlabel('Epoch [s]'); title('Altitude'); ylabel('[km]')

figure();
plot3(state(:, 1)./1E3, state(:, 2)./1E3, state(:, 3)./1E3, 'LineWidth', 2);
title('LRO orbit'); axis equal; grid on;

%% AUXILIARY FUNCTIONS
function [dx] = EOM_LRO_EPHEM_SRP_TEST(t, x, planetParams, C_mat, S_mat, A, M)
    %%                       EOM FUNCTION EPHEM
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 09/24/2025
    %
    %   Description: Compute the EoM using the Ephemerides model in
    %   cislunar space.
    %
    %   Input:
    %       t: time vector
    %       x: state vector [r1, r2, r3, v1, v2, v3, eta, bias vector]'
    %       planetParams: planet parameters 
    %                     [GM, Re, nmax, normalized]
    %       poleParams: pole parameters
    %                   [W, W0, RA, DEC]
    %       Cmat: SH C coefficients
    %       Smat: SM S coefficients
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%
    
    % number of states (pos + vel + bias)
    Nx = 6;
    
    % Planet constants
    GM_E = planetParams(1); GM_M = planetParams(2);
    R_E  = planetParams(3); R_M  = planetParams(4);
    
    % Get GM for the Sun [m^3/s^2]
    [GM_S] = cspice_bodvrd('SUN', 'GM', 1)*1E9; 
    [R_Sun] = cspice_bodvrd('SUN', 'RADII', 3).*1E3;      % Get R for the Sun  [m]

    % gravity computation params
    normalized = planetParams(6);
    n_max      = planetParams(5); 
    
    % compute Earth position. ref: J2000
    target = 'EARTH';
    et = t;                         % Convert UTC time to ephemeris time

    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';             % Set the observer to the MOON

    [Estate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Epos  = Estate(1:3)*1E3;                                            % [m]
    r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

    % compute Moon position. ref: MOON
    target = 'MOON';
    [Mstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Mpos  = Mstate(1:3)*1E3;                                            % [m]
    r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

    % compute Sun position. ref: MOON
    target = 'SUN';
    [Sstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Spos = Sstate(1:3)*1E3;                                             % [m]
    r3 = [x(1)-Spos(1);x(2)-Spos(2);x(3)-Spos(3)];                      % SC-Sun
    
    % relative vector from Earth & Sun to Moon
    r_ME = Mpos - Epos;
    r_MS = Mpos - Spos;

    % compute orientation
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et);

    % SRP acceleration
    F = shadow_function(R_Sun(1), R_M, Spos, r2);
    % F = 0;
    [aSRP, daSRP_dr, ~] = SRP(r3, 1, M, A);

    % compute gravity acceleration
    Cmat_E = C_mat{1};
    Smat_E = S_mat{1};
    [~, dU1, ddU1] = potentialGradient_nm(Cmat_E, Smat_E, 0, ...
                                                J2000_EARTH'*r1, R_E, GM_E, ...
                                                normalized);

    [~, dU1_T, ~] = potentialGradient_nm(Cmat_E, Smat_E, 0, ...
                                                J2000_EARTH'*r_ME, R_E, GM_E, ...
                                                normalized);

    Cmat_M = C_mat{2};
    Smat_M = S_mat{2};
    [~, dU2, ddU2] = potentialGradient_nm(Cmat_M, Smat_M, n_max, ...
                                                J2000_MOON'*r2, R_M, GM_M, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    dU1  = J2000_EARTH  * dU1;
    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';

    dU2  = J2000_MOON  * dU2;
    ddU2 = J2000_MOON  * ddU2  * J2000_MOON';
    
    % Sun acceleration on the S/C. Point mass
    dU3  = -GM_S * (r3 / (vecnorm(r3)^3));    

    % Tidial acceleration
    a_tidial_E = - J2000_EARTH  * dU1_T;
    a_tidial_S = GM_S * r_MS / (vecnorm(r_MS)^3);
    
    % total acceleration
    dU = dU2 + dU1 +  dU3 + a_tidial_E + a_tidial_S + F * aSRP;

    % compute gravity position partials
    T = ddU2 + ddU1 + F * daSRP_dr;
    
    % compute Jacobian
    J = compute_jacobian(T);

    % STM value 
    PHI = reshape(x(Nx+1:Nx+Nx*Nx), [Nx, Nx]);

    PHI_dot = J * PHI;  

    % differential equations
   dx =  [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3);
          reshape(PHI_dot, [Nx*Nx, 1])];
end

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