clear;
clc;
close all;
addpath('simplified_functions/');
addpath('data/');
format long g;
cspice_furnsh(...
    '/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')

%% INFORMATION MATRIX PER INCLINATION
% Compute the information matrix per inclination for a Lunar orbit
% Date: 10/10/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Body constants (Moon defaults; change as needed) ----
R   = 1737.4E3;                 % [m] 
mu  = 4.902800118E3;            % [km^3/s^2] 
[planetParams, Cmat_true, Smat_true] = load_universe();

% time
utc_start = '2025-05-20 00:00:00';
utc_stop  = '2025-05-20 02:00:00';
N         = 2000;           % number of samples
et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

% filter
fCut = 1E-4;
Fs = 1 / (et(2) - et(1));
d = designfilt('highpassiir', ...
    'FilterOrder', 12, ...
    'HalfPowerFrequency', fCut, ...
    'SampleRate', Fs);

% ---- Choose ONE set of orbital elements (radians) ----
a     = R/1E3 + 90;              % [km]
e_range = 0:1E-2:0.15;            % eccentricity [-]  
i_range = 0:1:90;                % inclination [deg]
Omega   = deg2rad(0);              % RAAN [rad]
omega   = deg2rad(0);              % arg. perigee [rad]
M       = deg2rad(0);              % mean anomaly [rad]

SNR = ones(6, length(i_range), length(e_range)) * NaN;

%% compute inclination = 0
tic;
incIdx = 1;
disp('Computing inclination ...  ' +string(i_range(incIdx)));
for eccIdx = 1:length(e_range)
    i = deg2rad(i_range(incIdx));
    e = e_range(eccIdx);
    
    % Compute initial conditions. X0
    [rECI, vECI] = coe2rv(a,e,i,Omega,omega,M,mu);
    X0 = [rECI;vECI].*1E3;  % [m]
    
    % integrate trajectory
    options = odeset('RelTol',1E-13,'AbsTol',1E-13);
    STM0 = reshape(eye(6,6), [36, 1]);
    [t, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
        Cmat_true, Smat_true), et, [X0; STM0], options);
    
    % compute RTN matrix
    BN_mat = ones(3*length(t), 3)*NaN;
    for k = 1:length(t)
        maxInd = 3 * k; minInd = maxInd - 2;
        [NB] = RTN2ECI(state(k, 1:3)', state(k, 4:6)');
        BN_mat(minInd:maxInd, :) = NB';
    
    end
    
    [SNR_dB, ~] = compute_SNR(t, state, Cmat_true, ...
    Smat_true, BN_mat, d);
    
    SNR(:, incIdx, eccIdx) = SNR_dB; 
end
dt = toc;
expected_time = dt * 90 / 60;
disp('      Expected computation time ... ' + string(expected_time) + ' min');

%% compute rest of coefficients
 for incIdx = 2:length(i_range)
     disp('Computing inclination ...  ' +string(i_range(incIdx)));
     for eccIdx = 1:length(e_range)
         i = deg2rad(i_range(incIdx)); e = e_range(eccIdx);
    
        % Compute initial conditions. X0
        [rECI, vECI] = coe2rv(a,e,i,Omega,omega,M,mu);
        X0 = [rECI;vECI].*1E3;  % [m]
        
        % integrate trajectory
        options = odeset('RelTol',1E-13,'AbsTol',1E-13);
        STM0 = reshape(eye(6,6), [36, 1]);
        [t, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
            Cmat_true, Smat_true), et, [X0; STM0], options);
        
        % compute RTN matrix
        BN_mat = ones(3*length(t), 3)*NaN;
        for k = 1:length(t)
            maxInd = 3 * k; minInd = maxInd - 2;
            [NB] = RTN2ECI(state(k, 1:3)', state(k, 4:6)');
            BN_mat(minInd:maxInd, :) = NB';
        end
    
        [SNR_dB, ~] = compute_SNR(t, state, Cmat_true, ...
        Smat_true, BN_mat, d);
    
        SNR(:, incIdx, eccIdx) = SNR_dB; 
     end
 end
 
figure()
tt = ["\Gamma_{RR}", "\Gamma_{RT}", "\Gamma_{RN}", "\Gamma_{TT}", ...
    "\Gamma_{TN}", "\Gamma_{NN}"];
[x, y] = meshgrid(e_range, i_range);
for k = 1:6
    subplot(2, 3, k)
    contourf(x, y, squeeze(SNR(k, :, :)))
    xlabel('Eccentrycity [-]')
    ylabel('inclination [deg]')
    colormap 'jet';
    colorbar();
    title(tt(k))
    % % clim([55 100])
end

%% FUNCTIONS
function [rI, vI] = coe2rv(a,e,i,Omega,omega,M,mu)
    % Elements -> ECI r,v (scalar)
    E = keplerE(M,e);
    f = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    p = a*(1 - e^2);
    r = p/(1 + e*cos(f));
    r_pf = [r*cos(f); r*sin(f); 0];
    v_pf = sqrt(mu/p)*[-sin(f); e+cos(f); 0];
    
    Q = rotz(Omega)*rotx(i)*rotz(omega);
    rI = Q*r_pf; vI = Q*v_pf;
end

function E = keplerE(M,e)
    % Solve Kepler's equation M = E - e sin E (scalar)
    M = mod(M, 2*pi);
    if e < 0.8, E = M; else, E = pi; end
    for k=1:50
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E  = E + dE;
        if abs(dE) < 1e-14, break; end
    end
end

function R = rotx(a)
    ca=cos(a); sa=sin(a);
    R = [1 0 0; 0 ca -sa; 0 sa ca];
    end

function R = rotz(a)
    ca=cos(a); sa=sin(a);
    R = [ca -sa 0; sa ca 0; 0 0 1];
end

function [SNR_dB, Y_filt] = compute_SNR(t, state, Cmat_true, ...
    Smat_true, BN_mat, d)

    % measurement weight
    sigma_noise = 1E-12;

    % Moon parameters
    [GM2] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    [Re2] = cspice_bodvrd('MOON', 'RADII', 3); % get Moon Radius [Km]
    GM2 = GM2*1E9;                             % [m^3 s^-2]
    Re2 = Re2*1E3;                             % [m]

    Y   = zeros(6, length(t));
    % compute information matrix
    for k = 1:length(t)
        % current position ECI
        r_ECI = state(k, 1:3)';

        % rotation matrix to RNT
        maxInd = 3 * k; minInd = maxInd - 2;
        BN     = BN_mat(minInd:maxInd, :);

        % rotation matrix
        frame_from = 'MOON_PA';
        frame_to   = 'J2000';
        J2000_MOON = cspice_pxform(frame_from, frame_to, t(k));

    
       % compute measurements at current time in RTN frame
       Cmat2 = Cmat_true{2};
       Smat2 = Smat_true{2};
       [~, ~, ddU_J2000] = potentialGradient_nm(Cmat2, Smat2, 160, ...
                                                J2000_MOON'*r_ECI, Re2(1), ...
                                                GM2, 1);
% %        [~, ~, ddU0_J2000] = potentialGradient_nm(Cmat2, Smat2, 0, ...
% %                                                 J2000_MOON'*r_ECI, Re2(1), ...
% %                                                 GM2, 1);
% %        ddU_J200 = J2000_MOON * (ddU_J2000 - ddU0_J2000) * J2000_MOON';
       ddU_J200 = J2000_MOON * (ddU_J2000) * J2000_MOON';

       ddU_RTN = BN * ddU_J200 * BN';

       % store measurements 
       Y(:, k) = [ddU_RTN(1,1); ddU_RTN(1,2); ddU_RTN(1,3); ddU_RTN(2,2);...
           ddU_RTN(2,3);ddU_RTN(3,3)];
    end

    % Filter measurements
    N      = length(Y(1, :));
    Y_filt = filtfilt(d, Y')';

    % compute SNR
    SNR = (1/ N) * sum(Y_filt.^2, 2)./(sigma_noise.^2);

    SNR_dB = 10 * log10(SNR);
end