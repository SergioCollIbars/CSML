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

%% ATTITUDE ERROR EFFECT IN 9:2 NRHO ORBIT
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

% Initial configuration
system = "EPHEM";
consider_cov = 0;
tmin = 0;                           % [rad]
tmax = 1*1.4968;                    % [rad]
frec = 1/1;                         % [Hz]
orientation = "RTN";           % Inertial / RTN / Earth / Sun
MC   = 2;

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME, DOM] = ...
    load_universe(system, [tmin, tmax], frec);

% load initial conditions
X0 = load_initCond(system, planetParams, TIME);

f_time = frec;                         % fixed integ. frec [Hz]
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);

% Orbit integrator
disp('Integrating trajectory ....')
options = odeset('RelTol',1E-13,'AbsTol',1E-13);
STM0 = reshape(eye(6,6), [36, 1]);
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0, 0), ...
    TIME, [X0; STM0], options);
disp('      DONE!');

% save integrated state
Xsc_ODE = state(:, 1)';
Ysc_ODE = state(:, 2)';
Zsc_ODE = state(:, 3)';

VXsc_ODE = state(:, 4)';
VYsc_ODE = state(:, 5)';
VZsc_ODE = state(:, 6)';

% compute orientation
BN_mat = nan(3 * length(TIME),  3);
for j = 1:length(TIME)
    et = TIME(j)./planetParams(3); % Convert UTC time to ephemeris time
    ref = 'J2000';
    abcorr = 'NONE';
    observer = '3';                 % Set the observer to the Earth-Moon barycenter
    
    % Moon position
    target = 'MOON';
    [stateM, ~] = cspice_spkezr(target, et, ref, abcorr, observer);
    posM = stateM(1:3)./planetParams(2) * 1E3;
    velM = stateM(4:6)./(planetParams(2)*planetParams(3)) * 1E3;

    target = 'EARTH';
    [stateM, ~] = cspice_spkezr(target, et, ref, abcorr, observer);
    posE = stateM(1:3)./planetParams(2) * 1E3;
    velE = stateM(4:6)./(planetParams(2)*planetParams(3)) * 1E3;
    
    target = 'SUN';
    [stateM, ~] = cspice_spkezr(target, et, ref, abcorr, observer);
    posS = stateM(1:3)./planetParams(2) * 1E3;
    velS = stateM(4:6)./(planetParams(2)*planetParams(3)) * 1E3;

    x = [Xsc_ODE(j);Ysc_ODE(j);Zsc_ODE(j);VXsc_ODE(j);VYsc_ODE(j);...
        VZsc_ODE(j)];
    p_sc_moon = x(1:3) - posM; % S/C - Moon
    v_sc_moon = x(4:6) - velM; % S/C - Moon

    J2000_RTN   = RTN2ECI(p_sc_moon, v_sc_moon);
    EARTH_J2000 = ECI2EARTH(posE,velE, x(1:3), x(4:6));
    SUN_J2000   = ECI2SUN(posS,x(1:3));
    EARTHSUN_J2000 = ECI2EARTHSUN(posE, posS, x(1:3));

    maxIndx = 3 * j; minIndx = maxIndx - 2;
    if(orientation == "RTN")
        BN_mat(minIndx:maxIndx, :) = J2000_RTN';
    elseif(orientation == "Inertial")
        BN_mat(minIndx:maxIndx, :) = eye(3);
    elseif(orientation == "Earth")
        BN_mat(minIndx:maxIndx, :) = EARTH_J2000;
    elseif(orientation == "Sun")
        BN_mat(minIndx:maxIndx, :) = EARTHSUN_J2000;
    end
end

% SC - MOON vector
et = TIME./planetParams(3); % Convert UTC time to ephemeris time
target = 'MOON';
[stateM, ~] = cspice_spkezr(target, et, ref, abcorr, observer);
posM = stateM(1:3, :)./planetParams(2) * 1E3;

r_SC_MOON = state(:, 1:3)' - posM;
figure()
plot(TIME, vecnorm(r_SC_MOON), 'LineWidth', 2); grid on;
title('S/C - MOON vector norm');


% compute nominal angular velocity
[angVel_vec] = angularVel_from_DCM(BN_mat, 1/frec);
[angAcc_vec] = compute_angAcc(angVel_vec, 1/frec);

% time to date converter
TIME = t';   % [-]
jd = 2451545 + TIME / planetParams(3) / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', ...
    'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';

date = humanReadableTime;

% plot norm of angular velocity
figure()
plot(date(2:end-2), vecnorm(angVel_vec(:, 2:end-2)), 'LineWidth', 2)
grid on; ylabel('[rad / s]');

% compute attitude error effects along NRHO
disp('Computing attitude error residuals ...')
arcsec          = 1 * sqrt(frec);                      % [arcseconds]
radians         = arcsec * pi / (180 * 3600);          % [rad]
% radians_per_sec = 0.002 * (pi/180) / sqrt(3600) / sqrt(1/frec);
radians_per_sec = 3E-8;

[b] = compute_FOMP(TIME./planetParams(3), radians_per_sec);
dt = 1/frec;                                           % seconds

deltaE_att = ones(6, length(date), MC) * NaN; deltaE_angVel = deltaE_att;
deltaE_angAcc = deltaE_att;
for mc = 1:MC
    disp('MC = ' + string(mc));
    Ath = normrnd(0, radians, [3,length(date)]);        % mean 0, std: radians
    Aw  = normrnd(0, radians_per_sec, [3,length(date)]); % mean 0, std: radians
    for j = 3:length(date)-2
        % current S/C position
        x = [Xsc_ODE(j);Ysc_ODE(j);Zsc_ODE(j);VXsc_ODE(j);VYsc_ODE(j);...
            VZsc_ODE(j)];
    
        % Cumpute gradiometer measurement
        [~, ~, ~, ~, ~, ~, ddU] = ...
            compute_sc_acceleration(TIME(j), x, planetParams, Cmat_true, Smat_true);
    
        % Attitude partials
        Y   = reshape(ddU, [9, 1]); 
        
        maxIndx = 3 * j; minIndx = maxIndx - 2;
        BN = BN_mat(minIndx:maxIndx, :);
    
        [Hrot] = compute_rotPartials_analy(Y, BN);
    
        % Measurement residual. Attitude error
        deltaE_att(:, j, mc) = ([Hrot(1:3, :);Hrot(5:6,:);Hrot(9, :)] * Ath(:, j))...
            .* (planetParams(3)^2); % [1/s^2]
    
        % Measurement residuals. Angular velocity error

        omega_dyad_nom     = dyad_operator(angVel_vec(:, j));
        omega_dyad_true    = dyad_operator(angVel_vec(:, j) + Aw(:, j));
        
        B_nom  = omega_dyad_nom  * omega_dyad_nom;
        B_true = omega_dyad_true * omega_dyad_true;
        B      = B_true - B_nom;
        deltaE_angVel(:, j, mc) = [B(1,1);B(1,2);B(1,3);...
            B(2,2);B(2,3);B(3,3)];

        omegaDot_dyad_nom   = dyad_operator(angAcc_vec(:, j));
        omega_dyad_true     = dyad_operator(angAcc_vec(:, j) + Aw(:, j)./dt);
        C = omega_dyad_true - omegaDot_dyad_nom;
        deltaE_angAcc(:, j, mc) = [C(1,1);C(1,2);C(1,3);...
            C(2,2);C(2,3);C(3,3)];
    end
end
disp('      DONE!')
disp('--------------------------------')

% load gravity errors
% % dY = load("disturbance_grav_error.mat").dY;

% plot error
figure()
idx = [1, 2, 3, 4, 5, 6];
lb = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", "\Gamma_{yy}", ...
    "\Gamma_{yz}", "\Gamma_{zz}"];
for k = 1:MC
    for j = 1:6
        subplot(2, 3, idx(j))

        if(k == 1)
            semilogy(date, ones(1, length(date)) * 3E-3 * sqrt(frec), ...
                'LineWidth', 2, 'Color','b'); hold all;
            semilogy(date, ones(1, length(date)) * 3E-3 * sqrt(frec), ...
                'LineWidth', 2, 'Color','r');
            semilogy(date, ones(1, length(date)) * 3E-3 * sqrt(frec), ...
                'LineWidth', 2, 'Color','g');
                        semilogy(date, ones(1, length(date)) * 3E-3 * sqrt(frec), ...
                'LineWidth', 2, 'Color','k');
        end

        semilogy(date, abs(deltaE_angVel(j, :, k))./1E-12, 'LineWidth', 2, 'Color', 'b', ...
            'LineStyle','none', 'Marker', '.', 'MarkerSize', 2); 
        semilogy(date, abs(deltaE_att(j, :, k))./1E-12, 'LineWidth', 2, 'Color', 'r', ...
            'LineStyle','none', 'Marker', '.', 'MarkerSize', 2);
% %         semilogy(date, abs(dY(j, :, k))./1E-9, 'LineWidth', 2, 'Color', 'g', ...
% %             'LineStyle','none', 'Marker', '.', 'MarkerSize', 2);
        if(k ==  1)
            xlabel('date')
            ylabel(lb(j) + '[mE]')
            grid on;
        end
    end
end
sgtitle('Observation error along NRHO orbit');
% % legend('angular velocity errors','orientation errors',...
% %     'gravity field errors','noise level');

% compute statistics
[AngVel_err_stats] = compute_stats(deltaE_angVel(:, 3:end-2, :)./1E-12);
[AngAcc_err_stats] = compute_stats(deltaE_angAcc(:, 3:end-2, :)./1E-12);
[Att_err_stats]    = compute_stats(deltaE_att(:, 3:end-2, :)./1E-12);

disp('Angular velocity errors stats')
display_stats(AngVel_err_stats);

disp('Angular acceleration errors stats')
display_stats(AngAcc_err_stats);

disp('Attitude errors stats')
display_stats(Att_err_stats);

figure(); titles= ["XX", "XY", "XZ", "YY", "YZ", "ZZ"];
colors = [
    0.0000  0.4470  0.7410   % blue
    0.8500  0.3250  0.0980   % orange
    0.9290  0.6940  0.1250   % yellow
    0.4940  0.1840  0.5560   % purple
    0.4660  0.6740  0.1880   % green
    0.3010  0.7450  0.9330   % cyan
    ];
for comp = 1:6
    subplot(6, 1, comp);
    for mc = 1:MC
        plot(date, squeeze(deltaE_att(comp, :, mc))./1E-12, ...
            'LineStyle', 'none', 'Marker','.', 'color', colors(comp, :));
        ylabel(titles(comp)); grid on; hold on;
    end
end
sgtitle('Atittude error induced residuals');
figure();
for comp = 1:6
    subplot(6, 1, comp);
    for mc = 1:MC
        plot(date, squeeze(deltaE_angVel(comp, :, mc))./1E-12, ...
            'LineStyle', 'none', 'Marker','.', 'color', colors(comp, :));
        ylabel(titles(comp)); grid on; hold on;
    end
end
sgtitle('Angular velocity error induced residuals');

figure();
Y_mag = squeeze(vecnorm(deltaE_att, 2, 1));
for k = 1:MC
    plot(date, Y_mag(:, k)./1E-12,...
        'Color', 'g', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10);
    hold on;
end
ylabel('milli-Eotvos')
grid on;
t_event = datetime(2020,1,15,09,09,32);   % example date/time
xline(t_event, '--', 'Color', 'r', 'LineWidth', 2);

t_start = datetime(2020,1,15,8,0,0);
t_end   = datetime(2020,1,15,10,0,0);
xlim([t_start t_end]);
title('Atittude errors impact on the gradiometer residuals');

figure();
Y_mag = squeeze(vecnorm(deltaE_angVel, 2, 1));
for k = 1:MC
    plot(date, Y_mag(:, k)./1E-12,...
        'Color', 'b', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 10);
    hold on;
end
ylabel('milli-Eotvos')
grid on;
t_event = datetime(2020,1,15,09,09,32);   % example date/time
xline(t_event, '--', 'Color', 'r', 'LineWidth', 2);

t_start = datetime(2020,1,15,8,0,0);
t_end   = datetime(2020,1,15,10,0,0);
xlim([t_start t_end]);
title('Angular velocity errors impact on the gradiometer residuals');

%% COMPUTE RMS AT PERILUNE PER COMPONENT
comp_lbl = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"];

% Attitude errors
disp('RMS PERILUNE. Attitude errors');
t_event = datetime(2020,1,15,09,09,32);   % perilune date
idx = date >= t_event - seconds(15*60) & date <= t_event + seconds(15*60);
for comp = 1:6
    data = squeeze(deltaE_att(comp, :, 1))./1E-12; % milli-Eotvos
    RMS_perilune = sqrt(mean(data(idx).^2, 'all'));
    disp(comp_lbl(comp) + ' ' + string(RMS_perilune))
end

% Angular velocity errors
disp('RMS PERILUNE. Angular velocity errors');
idx = date >= t_event - seconds(10*60) & date <= t_event + seconds(10*60);
for comp = 1:6
    data = squeeze(deltaE_angVel(comp, :, :))./1E-12; % milli-Eotvos
    RMS_perilune = sqrt(mean(data(idx,:).^2, 'all'));
    disp(comp_lbl(comp) + ' ' + string(RMS_perilune))
end

%% COMPUTE RMS AT APOLUNE PER COMPONENT

% Attitude errors
disp('RMS APOLUNE. Attitude errors');
t_event = datetime(2020, 1, 12, 04,27,00); % apolune date
idx = date >= t_event - seconds(60) & date <= t_event + seconds(60);
for comp = 1:6
    data = squeeze(deltaE_att(comp, :, :))./1E-12; % milli-Eotvos
    RMS_perilune = sqrt(mean(data(idx,:).^2, 'all'));
    disp(comp_lbl(comp) + ' ' + string(RMS_perilune))
end

% Angular velocity errors
disp('RMS APOLUNE. Angular velocity errors');
idx = date >= t_event - seconds(60) & date <= t_event + seconds(60);
for comp = 1:6
    data = squeeze(deltaE_angVel(comp, :, :))./1E-12; % milli-Eotvos
    RMS_perilune = sqrt(mean(data(idx,:).^2, 'all'));
    disp(comp_lbl(comp) + ' ' + string(RMS_perilune))
end

%% FUNCTIONS
function [] = display_stats(A)
    fprintf('      RMS[mE]        Max[mE]        Min[mE]\n');
    fprintf('-------------------------------------\n');
    
    for i = 1:size(A,1)
        fprintf('%8.2g   %8.2g   %8.2g\n', A(i,1), A(i,2), A(i,3));
    end
end
function [Y_stats] = compute_stats(Y)
   % RMS over time for each component and MC
    RMS_time_MC = sqrt(mean(Y.^2, 2, 'omitnan'));   % size: [6, 1, Mc]
    
    % Mean RMS across MC
    RMS_time_mean = mean(RMS_time_MC, 3);   % [6 x 1]

    Y_max = max(Y, [], [2 3]);   % [6 x 1]
    Y_min = min(Y, [], [2 3]);   % [6 x 1]

    
    Y_stats = [RMS_time_mean, Y_max, Y_min];
end

function [dUE, dUM, dUS, dUJ, dUEM, dUSRP, ddU] = compute_sc_acceleration(t, x, planetParams, C_mat, S_mat)
        [GM2] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
        [Re2] = cspice_bodvrd('MOON', 'RADII', 3); % get Sun Radius [Km]
        GM2 = GM2*1E9/(planetParams(2)^3 * planetParams(3)^2);
        Re2 = Re2*1E3./planetParams(2);

        [GM1] = cspice_bodvrd('EARTH', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
        [Re1] = cspice_bodvrd('EARTH', 'RADII', 3); % get Jupiter Radius [Km]
        GM1 = GM1*1E9/(planetParams(2)^3 * planetParams(3)^2);
        Re1 = Re1*1E3./planetParams(2);

        [GM3] = cspice_bodvrd('SUN', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
        [Re3] = cspice_bodvrd('SUN', 'RADII', 3); % get Sun Radius [Km]
        GM3 = GM3*1E9/(planetParams(2)^3 * planetParams(3)^2);
        Re3 = Re3*1E3./planetParams(2);

        [GM4] = cspice_bodvrd('5', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
        [Re4] = cspice_bodvrd('JUPITER', 'RADII', 3); % get Jupiter Radius [Km]
        GM4 = GM4*1E9/(planetParams(2)^3 * planetParams(3)^2);
        Re4 = Re4*1E3./planetParams(2);
        
        % gravity computation params
        n_max      = planetParams(6);
        normalized = planetParams(7);

        % compute Earth position. ref: J2000
        target = 'EARTH';
        et = t./planetParams(3);      % Convert UTC time to ephemeris time
        et2 = t./planetParams(3) + 1; % Convert UTC time to ephemeris time
        ref = 'J2000';
        abcorr = 'NONE';
        observer = '3';  % Set the observer to the Earth-Moon barycenter

        [Estate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Epos  = Estate(1:3)*1E3/planetParams(2);
        r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

        % compute Moon position. ref: J2000
        target = 'MOON';
        [Mstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Mpos  = Mstate(1:3)*1E3/planetParams(2);
        r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

        % compute Sun position. ref: J2000
        target = 'SUN';
        [Sstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Spos = Sstate(1:3)*1E3/planetParams(2);
        r3 = [x(1)-Spos(1);x(2)-Spos(2);x(3)-Spos(3)];                      % SC-Sun
        
        
        % EM barycenter acceleration. ref: J2000
        target = '3';
        [EMstate1, ~] = cspice_spkezr(target, et, ref, abcorr, '0');        % [Km & Km/s]
        [EMstate2, ~] = cspice_spkezr(target, et2, ref, abcorr, '0');       % [Km & Km/s]
        
        Svel  = EMstate1(4:6)*1E3;   % [m/s]
        Svel2 = EMstate2(4:6)*1E3;   % [m/s]
        
        At = (et2-et);
        Acc_EM = (Svel2- Svel)./At; % [m/s^2]
        Acc_EM = Acc_EM./(planetParams(2) * planetParams(3)^2); % [-]

        % compute Jupiter barycenter position. ref: J2000
        target = '5';
        [Jstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Jpos = Jstate(1:3)*1E3/planetParams(2);
        
        r4 = [x(1)-Jpos(1);x(2)-Jpos(2);x(3)-Jpos(3)];                      % SC-Sun

        % compute orientation
        frame_to   = 'J2000';
        frame_from = 'IAU_EARTH';
        J2000_EARTH = cspice_pxform(frame_from, frame_to, et);
        
        frame_from = 'MOON_ME';
        J2000_MOON = cspice_pxform(frame_from, frame_to, et);

        % compute gravity acceleration
        Cmat1 = C_mat{1};
        Smat1 = S_mat{1};
        [~, dU1, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                    J2000_EARTH'*r1, Re1(1), GM1, ...
                                                    normalized);
        Cmat2 = C_mat{2};
        Smat2 = S_mat{2};
        [~, dU2, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                    J2000_MOON'*r2, Re2(1), GM2, ...
                                                    normalized);

        Cmat3 = Cmat2;
        Smat3 = Smat2;
        Cmat3(2:end, :) = 0;
        Smat3(2:end, :) = 0;
        [~, dU3, ddUS] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                    r3, Re3(1), GM3, ...
                                                    normalized);

        [~, dU4, ddUJ] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                    r4, Re4(1), GM4, ...
                                                    normalized);

        % rotate back to inertial. Earth-Moon (EM) plane
        dU1  = J2000_EARTH * dU1;
        dU2  = J2000_MOON  * dU2;

        ddUE = J2000_EARTH * ddU1 * J2000_EARTH';
        ddUM = J2000_MOON  * ddU2 * J2000_MOON';

        % acceleration
        dUE = dU1;
        dUM = dU2;
        dUS = dU3;
        dUJ = dU4;
        dUEM = Acc_EM;

        % SRP force
        m = 1E3;                            % [Kg]
        A =  50/ (planetParams(2)^2);       % [-]
        eta  = 1.3;                        % [-]
        [dUSRP, ~, ~] = SRP(r3, eta, m, A, planetParams);

        % total gradiometer measurement
        ddU = ddUE + ddUM + ddUS + ddUJ;
end

function [M] = dyad_operator(x)
    % construct dyad operator
    M = [0, -x(3), x(2);...
        x(3), 0, -x(1);...
        -x(2), x(1), 0];
end


function [b] = compute_FOMP(time, sigma_b)
    N = length(time);
    b = zeros(3, N); 
    dt = time(2) - time(1); % [sec]
    tau = 5 * dt;
    var_b = sigma_b^2;
    qb = 2 * var_b / tau;
    Q = qb * tau / 2 * (1 - exp(-2*dt/tau));
    for j = 1:N-1
        w = normrnd(0, sqrt(Q), [3, 1]);
        b(:, j + 1)    = exp(-dt/tau) * b(:, j) + w;      % [1 / s^2]
    end
end