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
frec = 1/60;                        % [Hz]

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


% compute attitude error effects along NRHO
disp('Computing attitude error residuals ...')
arcsec          = 1;                                % [arcseconds]
radians         = arcsec * pi / (180 * 3600);       % [rad]
radias_per_sec  = 3E-7;                                 % [rad/s]
Ath = normrnd(0, radians, [3,length(date)]);        % mean 0, std: radians
Aw  = normrnd(0, radias_per_sec, [3,length(date)]); % mean 0, std: radians
[b] = compute_FOMP(TIME./planetParams(3), radias_per_sec);

deltaE_att = ones(6, length(date)) * NaN; deltaE_angVel = deltaE_att;
for j = 1:length(date)
    % current S/C position
    x = [Xsc_ODE(j);Ysc_ODE(j);Zsc_ODE(j);VXsc_ODE(j);VYsc_ODE(j);...
        VZsc_ODE(j)];

    % Cumpute gradiometer measurement
    [~, ~, ~, ~, ~, ~, ddU] = ...
        compute_sc_acceleration(TIME(j), x, planetParams, Cmat_true, Smat_true);

    % Attitude partials
    Y   = reshape(ddU, [9, 1]); BN = eye(3,3);
    [Hrot] = compute_rotPartials_analy(Y, BN);

    % Measurement residual. Attitude error
    deltaE_att(:, j) = ([Hrot(1:3, :);Hrot(5:6,:);Hrot(9, :)] * Ath(:, j))...
        .* (planetParams(3)^2); % [1/s^2]

    % Measurement residuals. Angular velocity error
    omega_dyad = dyad_operator(Aw(:, j));
% %     omega_dyad = dyad_operator(b(:, j));
    B = omega_dyad * omega_dyad;
    deltaE_angVel(:, j) = [B(1,1);B(1,2);B(1,3);B(2,2);B(2,3);B(3,3)];
end
disp('      DONE!')

% plot error
figure()
idx = [1, 2, 3, 4, 5, 6];
lb = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", "\Gamma_{yy}", ...
    "\Gamma_{yz}", "\Gamma_{zz}"];
for j = 1:6
    subplot(2, 3, idx(j))
    semilogy(date, ones(1, length(date)) * 1E-3, 'LineWidth', 2, 'Color','k')
    hold all;
    semilogy(date, abs(deltaE_angVel(j, :))./1E-9, 'LineWidth', 2, 'Color', 'b', ...
        'LineStyle','-')
    semilogy(date, abs(deltaE_att(j, :))./1E-9, 'LineWidth', 2, 'Color', 'r', ...
        'LineStyle','-')
    xlabel('date')
    ylabel(lb(j) + '[E]')
    grid on;
    ylim([1E-8, 1E-2]);
end
legend('mili-Eotvos', '\Delta Y');
sgtitle('Observation error along NRHO orbit');

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