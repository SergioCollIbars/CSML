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
%%                 QGG ODE INTEGRATOR VS SPICE TEST
% Description: test the ODE 113 integration accuracy vs the SPICE generated
% orbit. Tested orbit: Lunar gateway NRHO. 

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

% Initial configuration
system = "EPHEM";
consider_cov = 0;
tmin = 0;                           % [rad]
tmax = 6*1.4968;                    % [rad]
frec = 1/60;                        % [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME, DOM] = ...
    load_universe(system, [tmin, tmax], frec);

% load initial conditions
X0 = load_initCond(system, planetParams, TIME);

% numerical integration
tol = 1E-13;
for j = 1:length(tol)
    % current tolerance
    tolj = tol(j);

    % integrate trajectory
    options = odeset('RelTol',tolj,'AbsTol',tolj);
    STM0 = reshape(eye(6,6), [36, 1]);
    tStart = tic;
    [t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat_true, Smat_true, system, 0, {0,0}, 0), TIME, [X0; STM0], options);
    tEnd = toc(tStart);

    % save integrated state
    Xsc_ODE = state(:, 1)';
    Ysc_ODE = state(:, 2)';
    Zsc_ODE = state(:, 3)';

    VXsc_ODE = state(:, 4)';
    VYsc_ODE = state(:, 5)';
    VZsc_ODE = state(:, 6)';
end

TIME = t';   % [-]
jd = 2451545 + TIME / planetParams(3) / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', ...
    'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';

date = humanReadableTime;

% SPICE integration
time = TIME./planetParams(3);    % [s]
[sc_SPICE, ~] = cspice_spkezr('-60000', time, 'J2000', 'NONE', '3');
sc_SPICE(1:3, :) = sc_SPICE(1:3, :)*1E3/planetParams(2); 
sc_SPICE(4:6, :) = sc_SPICE(4:6, :)*1E3/(planetParams(2)*planetParams(3));

[MOON, ~] = cspice_spkezr('MOON', time, 'J2000', 'NONE', '3');
MOON(1:3, :) = MOON(1:3, :)*1E3/planetParams(2); 
MOON(4:6, :) = MOON(4:6, :)*1E3/(planetParams(2)*planetParams(3));

% motion of the primaries
[posE_ODE, posM_ODE] = compute_posPrimaries(TIME, planetParams, system);
cspice_kclear

% plot SPICE trajectory
figure()
plot3(sc_SPICE(1, :) , sc_SPICE(2, :) , sc_SPICE(3, :), 'LineWidth', 2)
hold on;
plot3(Xsc_ODE, Ysc_ODE, Zsc_ODE, 'LineWidth', 2)
axis equal;
grid on;

% plot trajectory components
figure()
subplot(2, 1, 1)
plot(date, vecnorm(sc_SPICE(1:3, :)), 'LineWidth', 2)
title('S/C ID : 6000')

subplot(2, 1, 2)
plot(date, vecnorm([Xsc_ODE;Ysc_ODE;Zsc_ODE]), 'LineWidth', 2)
title('Integrator')

% plot integration error
figure()
subplot(3, 2, 1)
semilogy(date, abs(sc_SPICE(1, :) - Xsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 2, 3)
semilogy(date, abs(sc_SPICE(2, :) - Ysc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 2, 5)
semilogy(date, abs(sc_SPICE(3, :) - Zsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')

subplot(3, 2, 2)
semilogy(date, abs(sc_SPICE(4, :) - VXsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vx [m/s]')
subplot(3, 2, 4)
semilogy(date, abs(sc_SPICE(5, :) - VYsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vy [m/s]')
subplot(3, 2, 6)
semilogy(date, abs(sc_SPICE(6, :) - VZsc_ODE)*planetParams(2)*planetParams(3), 'LineWidth', 2)
ylabel('Vz [m/s]')
xlabel('Epoch [s]')
sgtitle('Component error between SPCIE and filter integrator');


figure()
subplot(3, 1, 1)
plot(date, abs(MOON(1, :) - Xsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 1, 2)
plot(date, abs(MOON(2, :) - Ysc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 1, 3)
plot(date, abs(MOON(3, :) - Zsc_ODE)*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')
sgtitle('Position integrated traj. - Moon');


figure()
subplot(3, 1, 1)
plot(date, abs(sc_SPICE(1, :) - MOON(1,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('X [Km]')
subplot(3, 1, 2)
plot(date, abs(sc_SPICE(2, :) - MOON(2,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Y [Km]')
subplot(3, 1, 3)
plot(date, abs(sc_SPICE(3, :) - MOON(3,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('Z [Km]')
xlabel('Epoch [s]')
sgtitle('Position SPICE traj. - MOON');


figure()
plot(date, vecnorm(MOON(1:3, :) - [Xsc_ODE;Ysc_ODE;Zsc_ODE])*planetParams(2)/1E3, 'LineWidth', 2)
xlabel('Epoch [s]')
sgtitle('Position integrated traj. - Moon');

figure()
plot(date, vecnorm(sc_SPICE(1:3, :) - MOON(1:3,:))*planetParams(2)/1E3, 'LineWidth', 2)
ylabel('[Km]')
xlabel('date [days]')
sgtitle('Position norm SPICE traj. - MOON');


cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')
[sc_SPICE, ~] = cspice_spkezr('-60000', time, 'J2000', 'NONE', '3');
sc_SPICE(1:3, :) = sc_SPICE(1:3, :)*1E3/planetParams(2); 
sc_SPICE(4:6, :) = sc_SPICE(4:6, :)*1E3/(planetParams(2)*planetParams(3)); 
% compute spacecraft acceleration
dUE = ones(3, length(date)) * NaN; dUM = dUE; dUS = dUE; dUJ = dUE; 
dUEM = dUE; dUSRP = dUE;
for j = 1:length(date)
    x = sc_SPICE(1:6, j)';
    [dUE(:, j), dUM(:, j), dUS(:, j), dUJ(:, j), dUEM(:, j), dUSRP(:, j), ~] = ...
        compute_sc_acceleration(time(j)*planetParams(3), x, planetParams,...
        Cmat_true, Smat_true);
end

dU = dUE + dUM + dUS + dUJ + dUSRP - dUEM;

figure()
scale = planetParams(2) * planetParams(3)^2;
semilogy(date, vecnorm(dU)*scale, date, vecnorm(dUE)*scale, date, ...
    vecnorm(dUM)*scale, date, vecnorm(dUS)*scale, date, vecnorm(dUJ)*scale,...
    date, vecnorm(dUSRP)*scale, 'LineWidth', 2)
xlabel('Time [days]')
ylabel('[m/s^2]')
title('Acceleration in NRHO orbit. Earth-Moon barycenter frame')
legend('total', 'Earth', 'Moon', 'Sun', 'Jupiter', 'SRP')

% compute attitude error effects along NRHO
arcsec = 1;
Ath = arcsec * pi / (180 * 3600).*[1;1;1];
deltaE = ones(6, length(date)) * NaN;
deltaX = ones(1, length(date)) * NaN;
signal = deltaE;
for j = 1:length(date)
    x = [Xsc_ODE(j);Ysc_ODE(j);Zsc_ODE(j);VXsc_ODE(j);VYsc_ODE(j);...
        VZsc_ODE(j)];
    [Hrot, ddU] = compute_rotErr(time(j)*planetParams(3), x, planetParams, ...
        Cmat_true, Smat_true);
    [hpos] = compute_posErr(time(j)*planetParams(3), x(1:3), planetParams, ...
        Cmat_true, Smat_true);
    Hpos = [hpos(1:3, :);hpos(5:6,:);hpos(9, :)];
    deltaE(:, j) = ([Hrot(1:3, :);Hrot(5:6,:);Hrot(9, :)] * Ath).* (planetParams(3)^2); % [1/s^2]
    dY = deltaE(:, j)./ (planetParams(3)^2);                                            % [-]
    % % dY = 1E-3 * ones(6, 1) * 1E-9 / (planetParams(3)^2);
    deltaX(j) = vecnorm(inv(Hpos'*Hpos) * (Hpos' * dY)) * planetParams(2);              % [m]   
    signal(:, j) = [ddU(1:3, :);ddU(5:6,:);ddU(9, :)].* (planetParams(3)^2);            % [1/s^2]
end

figure()
semilogy(date, vecnorm(deltaE)./1E-9, date, deltaX)
xlabel('date [days]')
legend('signal error [E]', 'position error [m]')

figure()
idx = [1, 3, 5];
for j = 1:3
    subplot(2, 3, idx(j))
    semilogy(date, abs(signal(j, :))./1E-9, 'LineWidth', 2, 'Color','k')
    hold on;
    semilogy(date, abs(deltaE(j, :))./1E-9, 'LineWidth', 2, 'Color', 'r', 'LineStyle','-')
    xlabel('date [days]')
    ylabel('[Eotvos]')
end

idx = [1, 3, 5, 2, 4, 6];
for j = 4:6
    subplot(2, 3, idx(j))
    semilogy(date, abs(signal(j, :))./1E-9, 'LineWidth', 2, 'Color','k')
    hold on;
    semilogy(date, abs(deltaE(j, :))./1E-9, 'LineWidth', 2, 'Color', 'r', 'LineStyle','-')
    xlabel('Epoch [s]')
    ylabel('[Eotvos]')
end
legend('Signal', '\Delta E')
sgtitle('Observation error compared to gradiometer signal along NRHO orbit')

figure()
idx = [1, 2, 3, 4, 5, 6];
lb = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", "\Gamma_{yy}", ...
    "\Gamma_{yz}", "\Gamma_{zz}"];
for j = 1:3
    subplot(2, 3, idx(j))
    semilogy(date, ones(1, length(date)) * 1E-3, 'LineWidth', 2, 'Color','k')
    hold on;
    semilogy(date, abs(deltaE(j, :))./1E-9, 'LineWidth', 2, 'Color', 'r', 'LineStyle','-')
    xlabel('date')
    ylabel(lb(j) + '[E]')
    grid on;
end

Nt = length(date);
for j = 4:6
    subplot(2, 3, idx(j))
    semilogy(date, ones(1, length(date)) * 1E-3, 'LineWidth', 2, 'Color','k')
    hold on;
    semilogy(date, abs(deltaE(j, :))./1E-9, 'LineWidth', 2, 'Color', 'r', 'LineStyle','-')
    xlabel('date')
    ylabel(lb(j) + ' [E]')
    grid on;
end
legend('mili-Eotvos', '\Delta Y');
sgtitle('Observation error along NRHO orbit');

% close SPICE
cspice_kclear
%% FUNCTIONS

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

function [Hrot, T] = compute_rotErr(t, x, planetParams, C_mat, S_mat)
    % define rotational partials
    Hrot = ones(9, 3) * NaN;

    % compute acceleration at current epoch
    [~, ~, ~, ~, ~, ~, ddU] = compute_sc_acceleration(t, x, ...
        planetParams, C_mat, S_mat);
    eps = 1E-6;
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;

        Atpos = At./2;
        Atneg = - At./2; 
 
        % Rotation matrix
        [Rpos] = rotationMatrix(Atpos(1), Atpos(2), Atpos(3), [3, 2, 1]);
        [Rneg] = rotationMatrix(Atneg(1), Atneg(2), Atneg(3), [3, 2, 1]);

        ddUpos = Rpos' * ddU * Rpos;
        ddUneg = Rneg' * ddU * Rneg;
        H = (ddUpos - ddUneg)./(vecnorm(Atpos-Atneg));
        
      Hrot(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
        H(3,1); H(3,2);H(3,3)];
    end

    T = [ddU(1,1);ddU(1,2);ddU(1,3);ddU(2,1);ddU(2,2);ddU(2,3);...
        ddU(3,1); ddU(3,2);ddU(3,3)];
end

function [Hpos] = compute_posErr(t, x, planetParams, C_mat, S_mat)
    % define rotational partials
    Hpos = ones(9, 3) * NaN;

    eps = 1E-6;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        Arpos = x + Ar./2;
        Arneg = x - Ar./2; 

        % compute acceleration at current epoch
        [~, ~, ~, ~, ~, ~, ddU_pos] = compute_sc_acceleration(t, Arpos, ...
        planetParams, C_mat, S_mat);

        [~, ~, ~, ~, ~, ~, ddU_neg] = compute_sc_acceleration(t, Arneg, ...
        planetParams, C_mat, S_mat);
 
        H = (ddU_pos - ddU_neg)./(vecnorm(Arpos-Arneg));
        
      Hpos(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
        H(3,1); H(3,2);H(3,3)];
    end
end
