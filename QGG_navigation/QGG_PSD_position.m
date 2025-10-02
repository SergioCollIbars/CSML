clear; 
clc;
close all;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

% Initial configuration
system = "CR3BP";    % options: 2BP, CR3BP, F2BP, FCR3BP, EPHEM
frame  = "inertial"; % options: inertial, autonomous, RTN

% time parameters
T_orb = 1.3817;                        % [rad]
tmin = 0;                              % [rad]
tmax = 9 * T_orb;                      % [rad]
frec = 1/30;                            % [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME] = ...
    load_universe(system, [tmin, tmax], frec);

f_time = frec;                         % fixed integ. frec [Hz]
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);

% load initial conditions. inertial frame (baricenter centered)
X0 = load_initCond(system, planetParams, TIME);

%% integrate trajectory
disp('Computing trajectory ...')
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
STM0 = reshape(eye(6,6), [36, 1]);

[t, state_inertial] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0, 0), TIME, ...
    [X0; STM0], options);
TIME = t;

% select frame rotation
if(frame == "inertial")
    pos = state_inertial;
    I = eye(3); N =length(TIME); BN_mat = ones(3 * N, 3) * NaN;
    for k = 1:N, maxInd = 3*k; minInd = maxInd - 2; BN_mat(minInd:maxInd, :) = I; end
elseif(frame == "autonomous")
    [state] = rotate_2_autonomousFrame(state_inertial(:, 1:3)', TIME);
    pos = state';
    I = eye(3); N =length(TIME); BN_mat = ones(3 * N, 3) * NaN;
    for k = 1:N, maxInd = 3*k; minInd = maxInd - 2; BN_mat(minInd:maxInd, :) = I; end
elseif(frame == "RTN")
    [r_ECI, v_ECI] = compute_Moon_relCoord(planetParams(1), ...
        state_inertial, TIME);
    [~, pos_RTN, BN_mat] = rotate_2_RTNFrame([], r_ECI, v_ECI, TIME);
    pos = pos_RTN';
end

% Plot orbit
plot_orbit(pos, TIME./planetParams(3), planetParams)

%% Select gravity content
% % planetParams(6) = 2;
% % 
% % A = Cmat_true{1,1}; B = Cmat_true{1,2};
% % C = Smat_true{1,1}; D = Smat_true{1,2}; 
% % for k = 1:planetParams(6)
% %     A(k, :) = A(k, :).*0; B(k, :) = B(k, :).*0;
% %     C(k, :) = C(k, :).*0; D(k, :) = D(k, :).*0;
% % end
% % Cmat_true{1,1} = A; Cmat_true{1,2} = B;
% % Smat_true{1,1} = C; Smat_true{1,2} = D;

%% compute gradiometer signal (Inertial frame [1/s^2])
disp('Computing gravity gradient ... ')
[T] = computeGrad_signal(state_inertial, planetParams, poleParams,...
    TIME, Cmat_true, Smat_true);
% select frame rotation
if(frame == "autonomous")
    [T2] = rotate_2_autonomousFrame(T, TIME);
    T = T2;
elseif(frame == "RTN")
    [T2, ~] = rotate_2_RTNFrame(T, r_ECI, v_ECI, TIME);
    T = T2;
end

%% compute attitude partials
disp('Computing rotation partials ... ')
[Hrot_mat] = compute_rot_partial_matrix(T, TIME, BN_mat);

%% compute position partials
disp('Computing position partials ... ')
[Hpos_mat] = compute_pos_partial_matrix(planetParams, poleParams, ...
    state_inertial, Cmat_true, Smat_true, TIME, BN_mat);

% plot gravity gradients
plot_gravityGradients(T, TIME./(planetParams(3)*86400))

%% compute and plot PSD for the gravity signal
disp('Computing PSD of the signals ... ')
[PSD_T, f_T] = compute_signal_PSD(frec, T_orb./planetParams(3), T./1E-9);
signal = [PSD_T(1:3, :); PSD_T(5:6, :); PSD_T(9, :)];
plot_PSD_values(f_T, signal, "Gravity tensor signal");
PSD_T = signal;

% % %% compute and plot PSD for the gravity position sensitivity 
% % % % % % % % % % % % (x direction error)
% % H_x = squeeze(Hpos_mat(:, 1, :)./1E-9);
% % [PSD_Hx, f_H] = compute_signal_PSD(frec, T_orb./planetParams(3),  H_x);
% % plot_PSD_values(f_H, PSD_Hx, "X direction sensitivity");
% % % % % % % % % % % % (y direction error)
% % H_y = squeeze(Hpos_mat(:, 2, :)./1E-9);
% % [PSD_Hy, f_H] = compute_signal_PSD(frec,T_orb./planetParams(3), H_y);
% % plot_PSD_values(f_H, PSD_Hy, "Y direction sensitivity");
% % % % % % % % % % % % (z direction error)
% % H_z = squeeze(Hpos_mat(:, 3, :)./1E-9);
% % [PSD_Hz, f_H] = compute_signal_PSD(frec,T_orb./planetParams(3), H_z);
% % plot_PSD_values(f_H, PSD_Hz, "Z direction sensitivity");

%% compute and plot signal Spectrogram
disp('Computing Spectrogram of the signals ... ')
[F, Tt, SPCT] = compute_spectrogram(frec, T_orb./planetParams(3), T./1E-9);
vols = {SPCT(:, :, 1:3), SPCT(:, :, 5:6), SPCT(:, :, 9)};   
signal = cat(3, vols{:});                               
plot_spectrogram(F, Tt,  signal, T_orb./planetParams(3), ...
    "Spectrogram gravity tensor");

%% compute Fisher Information assuming white noise
disp('Computing Fisher Information FFT ... ')
Sn = (1E-3)^2 * ones(1, 6); % [E^2 / Hz]
[f, I, ~] = compute_FFT_info(Hpos_mat./1E-9, frec, Sn);
pName = {'X','Y','Z'}; tt = "Fisher information, position state";
[~] = plot_Information(I, f, tt, pName);

[f, I, ~] = compute_FFT_info(Hrot_mat./1E-9, frec, Sn);
pName = {'Yaw','Pithc','Roll'}; tt = "Fisher information, attitude state";
[~] = plot_Information(I, f, tt, pName);

% clear kernels
cspice_kclear
disp('DONE!')

%% FUNCTIONS

function [T] = computeGrad_signal(state, planetParams, poleParams, ...
    TIME, C_mat, S_mat)
    T = zeros(9, length(TIME));

    % S/C position
    pos = state(:, 1:3)'.*planetParams(2);      % [m]

    % Compute gravity gradient [1/s^2]
    for k = 1:length(TIME)
        % compute position Earth and Moon
        [posE, posM] = computePos_circular(TIME(k), planetParams);

        [ddU] = compute_gradiometer_measurements(pos(:, k), posE, posM, ...
                    planetParams, poleParams, C_mat, S_mat, TIME(k));
        T(:, k) = reshape(ddU', [9,1]);
    end
end

function [pos_earth, pos_moon] = computePos_circular(t, planetParams)
    mu = planetParams(1);
    M  = t;

    % Earth position
    pos_earth = -[mu*cos(M);mu*sin(M);0];           % [-]

    % Moon position
    pos_moon  = -[(mu-1)*cos(M);(mu-1)*sin(M); 0];  % [-]

    % bacl to meters
    pos_earth = pos_earth * planetParams(2);        % [m]
    pos_moon  = pos_moon  * planetParams(2);        % [m]
end

function [ddU] = compute_gradiometer_measurements(state, posE, posM, ...
    planetParams, poleParams, C_mat, S_mat, t)
    % rotation from Earth-Moon planet to J2000
    i_EM = deg2rad(0);   % [rad]
    EM_N = rotationMatrix(0, 0, i_EM, [1, 1, 1]);

    % extract rotation parameters
    RA_E = poleParams(1);            % RA Earth [rad]
    DEC_E = poleParams(2);           % DEC Earth [rad]
    W0_E = poleParams(3);            % prime meridian Earth [rad]
    WDot_E = poleParams(4);          % ang. velocity Earth [rad/s]

    RA_M = poleParams(5);            % RA Moon [rad]
    DEC_M = poleParams(6);           % DEC Moon [rad]
    W0_M = poleParams(7);            % prime meridian Moon [rad]
    WDot_M = poleParams(8);          % ang. velocity Moon [rad/s]

    Wt_E = WDot_E * t / planetParams(3) + W0_E;
    Wt_M = WDot_M * t / planetParams(3) + W0_M;
    ACAF1_N = rotationMatrix(pi/2 + RA_E, pi/2 - DEC_E, Wt_E, [3, 1, 3]);
    ACAF2_N = rotationMatrix(pi/2 + RA_M, pi/2 - DEC_M, Wt_M, [3, 1, 3]);
    
    ACAF1_EM = ACAF1_N * EM_N';
    ACAF2_EM = ACAF2_N * EM_N';

    % extract planet paramters (non-dimensional units)
    GM_E = planetParams(8);     % [m^3 s^-2]
    GM_M = planetParams(9);     % [m^3 s^-2]
    Re_E = planetParams(4);     % [m]
    Re_M = planetParams(5);     % [m]
    
    n_max      = planetParams(6);
    normalized = planetParams(7);

    relE = state(1:3) - posE;
    relM = state(1:3) - posM;
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, ~, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                ACAF1_EM*relE, Re_E, GM_E, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, ~, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                ACAF2_EM*relM, Re_M, GM_M, ...
                                                normalized);
    
    % Gravity gradient
    ddU =  (ACAF1_EM' * ddU1 * ACAF1_EM) + (ACAF2_EM' * ddU2 * ACAF2_EM);
end

function [] = plot_gravityGradients(T, time)
    % plotting options
    cl = 'g';
    lw = 2;
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];
    signal = [T(1:3, :); T(5:6, :); T(9, :)]./1E-9; % [Eotvos]

    figure()
    for k = 1:6
        subplot(2, 3, k)
        plot(time, signal(k, :), 'LineWidth', lw, 'Color', cl)
        xlabel('Time [days]')
        ylabel('Eotvos')
        title(tt(k));
        grid on;
    end
end

function [H] = compute_pos_partials_CR3BP(planetParams, poleParams, x,...
    C_mat, S_mat, posE, posM, t, BN)
    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [ACI]
        rneg = x - Ar./2;   % [ACI]

        [ddU_pos] = compute_gradiometer_measurements(rpos, posE, posM, ...
                        planetParams, poleParams, C_mat, S_mat, t);

        [ddU_neg] = compute_gradiometer_measurements(rneg, posE, posM, ...
                        planetParams, poleParams, C_mat, S_mat, t);

        Ht = ((BN * ddU_pos * BN') - (BN * ddU_neg * BN'))./...
            (vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end

function [Hpos_mat] = compute_pos_partial_matrix(planetParams, poleParams, ...
    state, C_mat, S_mat, TIME, BN_mat)
    Hpos_mat = ones(6,3,length(TIME)) * NaN;

     % S/C position
    pos = state(:, 1:3)'.*planetParams(2);      % [m]

    for k = 1:length(TIME)
         maxInd = 3* k; minInd = maxInd - 2;
         BN = BN_mat(minInd:maxInd, :); 
         % compute position Earth and Moon
        [posE, posM] = computePos_circular(TIME(k), planetParams);

        [H] = compute_pos_partials_CR3BP(planetParams, poleParams, pos(:, k),...
            C_mat, S_mat, posE, posM, TIME(k), BN);
        Hpos_mat(:, :, k) = H;
    end
end

function [Hrot_mat] = compute_rot_partial_matrix(T, TIME, BN_mat)
    Hrot_mat = ones(6,3,length(TIME)) * NaN;

    for k = 1:length(TIME)
        maxInd = 3 * k; minInd = maxInd - 2;
        BN  = BN_mat(minInd:maxInd, :);

        [hrot] = compute_rotPartials_analy(T(:, k), BN);
         Hrot_mat(:, :, k) = [hrot(1, :);hrot(2, :);hrot(3, :);hrot(5, :);...
            hrot(6,:);hrot(9, :)];
    end
end

function [PSD, frec] = compute_signal_PSD(fs, T_orbit, signal)
    % compute the PSD
    Lseg = max(8, round(T_orbit * fs));   
    window = hann(Lseg, "periodic");
    nWindow = length(window);
    noverlap = floor(0.5 * Lseg);  % 75%
    nfft = Lseg;
    
    % number of signal chanels
    Nc = length(signal(:, 1));

    % PSD output
    PSD = ones(Nc, nWindow/2 + 1);
    frec = PSD;
    for k =1:Nc
        % Compute PSD
        [PSD(k, :), frec(k, :)] = ...
            pwelch(signal(k, :), window, noverlap, nfft, fs, "onesided");
    end
end

function [F, T, PSD] = compute_spectrogram(fs, T_orbit, signal)
    % number of signal chanels
    Nc = length(signal(:, 1));
    
    % Spectrogram parameters
    L_orb   = round(T_orbit*fs);
    frac = 1/8;                               % try 1/4–1/8 orbit windows
    Lseg = max(128, round(frac*L_orb));
    window  = hann(Lseg,"periodic");
    hop  = max(1, round(0.05*Lseg));           % 95% overlap
    noverlap = Lseg - hop;
    nfft  = max(Lseg, 4*Lseg);
    
     [~,~,t,~] = spectrogram(signal(1, :), window, noverlap,...
            nfft, fs, "yaxis");  % PSD units: unit^2/Hz
     Nt  =length(t);

    % PSD output
    PSD = ones(Lseg*2 + 1, Nt, Nc) * NaN;
    F   = ones(Nc, Lseg*2 + 1) * NaN;
    T   = ones(Nc, Nt) * NaN;
    for k =1:Nc
        % compute Spectrogram
        [~,F(k, :),T(k, :),PSD(:, :, k)] = spectrogram(signal(k, :), window, noverlap,...
            nfft, fs, "yaxis");  % PSD units: unit^2/Hz
    end
end

function [] = plot_PSD_values(frec, PSD, sgtt)
    % plotting options
    cl = 'r';
    lw = 2;
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];

    % Number of chanels
    Nc = length(PSD(:, 1));

    figure()
    for k = 1:Nc
        subplot(2, 3, k)
        loglog(frec, PSD(k, :), 'LineWidth', lw, 'Color', cl)
        xlabel('Frequency [HZ]')
        ylabel('E^2 / Hz')
        % % ylabel('E^2 /m^2 / Hz')
        title(tt(k));
        grid on;
    end
    sgtitle(sgtt);
end

function [] = plot_orbit(state, time, planetParams)
    figure()
    tt = ["x", "y", "z"];
    pos = state(:, 1:3).*planetParams(2)./1E3;  % [km]
    for k = 1:3
        subplot(3, 1, k)
        plot(time./86400, pos(:, k)', 'LineWidth', 2)
        xlabel('Time [days]')
        ylabel('[Km]')
        title(tt(k))
        grid on;
    end
    
    figure()
    plot(time./86400, vecnorm(pos(:, 1:3)'), 'LineWidth', 2)
    xlabel('Time [days]')
    ylabel('[Km]')
    title('Orbit radius norm')

    figure()
    plot3(state(:, 1), state(:, 2), state(:, 3), 'LineWidth', 2)
    grid on;
    axis equal;
    title('Orbit in Inertial frame');
end

function [] = plot_spectrogram(F, T, PSD, T_orbit, sgtt)
    % Signal chanels
    Nc = length(PSD(1, 1, :));
    tt = ["xx", "xy", "xz", "yy", "yz", "zz"];
    figure();
    for k = 1:Nc
        subplot(2, 3, k)
        x = squeeze(PSD(:, :, k));
        imagesc(T(k, :)./T_orbit, F(k, :), 10*log10(x)); axis xy; ...
            colormap turbo; colorbar;
        xlabel('Orbit Number'); ylabel('Frequency [Hz]');
        title(tt(k))
    end
    sgtitle(sgtt)

    figure()
    x = squeeze(PSD(:, :, 1)); % [E^2 / Hz]
    pcolor(T(1, :)./T_orbit, F(1, :), 10*log10(x)); 
    axis xy; colormap turbo; colorbar;
    shading flat
    set(gca, 'YScale', 'log');    % Make y-axis log
    xlabel('Orbit fraction'); ylabel('Frequency [Hz]');
    title('Gradiometer spectrogram in xx direction')

    clim([-60, 40])
    ylim([0, 8E-3]);
end

function [var] = rotate_2_autonomousFrame(x, TIME)
    Nt = length(TIME);
    if(length(x(:, 1)) == 3)    % Position vector rotation
        var = ones(3, Nt) * NaN;
        for k = 1:Nt
            NB = [cos(TIME(k)), -sin(TIME(k)), 0;...
                  sin(TIME(k)), cos(TIME(k)), 0;...
                   0, 0, 1];
    
            var(:, k) = NB' * x(:, k);
        end
    elseif(length(x(:, 1)) == 9) % Gravity tensor rotation
        var = ones(9, Nt) * NaN;
        for k = 1:Nt
            NB = [cos(TIME(k)), -sin(TIME(k)), 0;...
                  sin(TIME(k)), cos(TIME(k)), 0;...
                   0, 0, 1];
    
            var(:, k) = reshape(NB' * reshape(x(:, k), [3,3]) * NB, [9,1]);
        end
    end
end

function [var, pos_RTN, BN] = rotate_2_RTNFrame(x, r_ECI, v_ECI, TIME)
    Nt = length(TIME);
    pos_RTN = ones(3, Nt) * NaN;
    var     = ones(9, Nt) * NaN; 
    BN      = ones(3*Nt, 3) * NaN;
    for k = 1:Nt
        r = r_ECI(:, k); v = v_ECI(:, k);
        n = cross(r, v);
    
        R = r / vecnorm(r);
        N = n / vecnorm(n);
        T = cross(N, R);
    
        RTN_ECI = [R'; T'; N'];

        maxInd = 3 * k; minInd = maxInd -2;
        BN(minInd:maxInd, :) = RTN_ECI;

        pos_RTN(:, k) = RTN_ECI * r;
        if(~isempty(x))
            var(:, k) = reshape(RTN_ECI * reshape(x(:, k), [3,3]) ...
                * RTN_ECI', [9,1]);
        end
    end
end

function [r_rel_ECI, v_rel_ECI] = compute_Moon_relCoord(mu, ...
    state_inertial, TIME)
    Nt  =length(TIME);
    r_ECI = state_inertial(:, 1:3)'; v_ECI = state_inertial(:, 4:6)';
    r_rel_ECI = r_ECI.*NaN; v_rel_ECI = v_ECI.*NaN;
    for k = 1:Nt
        % relative position w.r.t the Moon
        rM = -[(mu-1);0; 0]; 
        vM = [0;0;0];                      % Autonomous frame coordinates, [-]

        [rM_ECI, vM_ECI] = rotate2inertial(rM, vM, TIME(k), 1);  % [-] and [-]

        r_rel_ECI(:, k) = r_ECI(:, k) - rM_ECI; ...
        v_rel_ECI(:, k) = v_ECI(:, k) - vM_ECI;
    end
end

function [f, Iprime, Hf] = compute_FFT_info(Ht, fs, Sn)
    % Ht : 6 x 3 x Nt   (channels x params x time)
    % fs : sampling rate [Hz]
    % Sn : noise CSD. Any of:
    %      - Nf x 6              (diagonal, freq-dependent; one-sided)
    %      - 1  x 6  or 6 x 1    (diagonal, white)
    %      - 6  x 6              (full, white)
    %      - 6  x 6  x Nf        (full, freq-dependent; one-sided)
    %
    % Returns:
    %   f      : Nf x 1  frequency grid [Hz] (one-sided, f >= 0)
    %   Iprime : Nf x 3 x 3  Fisher information **per frequency bin** (no integration)
    %   Hf     : 6 x 3 x Nf  Fourier transform of Ht (continuous-time scaling)

    [ny, nth, Nt] = size(Ht);     % ny=6, nth=3
    dt = 1/fs;
    df = fs / Nt;

    % Continuous-time FT along time (dim 3)
    Hf_full = dt * fft(Ht, [], 3);     % 6 x 3 x Nt

    % Keep one-sided spectrum (assumes real Ht over time)
    Nf = floor(Nt/2) + 1;
    Hf = Hf_full(:,:,1:Nf);            % 6 x 3 x Nf
    f  = (0:Nf-1).' * df;              % Hz

    % Build per-frequency Fisher: I'(f_k) = H^*(f_k) S_n^{-1}(f_k) H(f_k)
    Iprime = zeros(Nf, nth, nth);
    for k = 1:Nf
        Hk = Hf(:,:,k);                % 6 x 3

        % Inverse noise weighting
        if nargin < 3 || isempty(Sn)
            W = eye(ny);
        else
            sz = size(Sn);
            if isequal(sz, [Nf, ny])           % diagonal, freq-dependent
                W = diag(1 ./ Sn(k,:));
            elseif (isequal(sz, [1, ny]) || isequal(sz, [ny, 1])) % diagonal, white
                W = diag(1 ./ Sn(:).');
            elseif isequal(sz, [ny, ny])       % full, white
                W = inv(Sn);
            elseif isequal(sz, [ny, ny, Nf])   % full, freq-dependent
                W = inv(Sn(:,:,k));
            else
                error('Sn has unexpected size.');
            end
        end

        Iprime(k,:,:) = Hk' * W * Hk;          % 3 x 3
    end
end

function [Idiag] = plot_Information(Iprime, f, tt, PName)
    Nf  = length(Iprime(:, 1, 1));
    nth = length(Iprime(1, :, 1));
    % Extract the diagonal (one curve per parameter)
    Idiag = zeros(Nf, nth);
    for k = 1:Nf
        Idiag(k,:) = real(diag(squeeze(Iprime(k,:,:))));
    end
    
    % ---- Plot (log axes recommended) ----
    paramNames = PName;             % adjust names if needed
    figure; 
    loglog(f, Idiag(:,1), f, Idiag(:,2), f, Idiag(:,3), 'LineWidth', 2)
    grid on
    xlabel('Frequency [Hz]')
    ylabel('1/param^2 per Hz')
    legend(paramNames, 'Location','best')
    title(tt)
end
