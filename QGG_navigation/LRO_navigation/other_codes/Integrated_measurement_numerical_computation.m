clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);
%% COMPUTE GRADIOMETER INTEGRATED MEASUREMENT NUMERICALLY
% Description: Using 2BP dynamics, compute the integrated grav. gradiometer
%  measurement and compare different numerical methods to approximate it.
%
% Author: Sergio Coll-Ibars
% Date: 06/05/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Planet constants
GM = 4902.800118 * 1E9;   % Moon [m^3/s^2]
R  = 1737.4*1E3;          % m

%% Initial condition given in orbital elements (oe)
% semi-major axis [m]
a = R + 2E4;

% eccentricity [-]
e = 0.8;

% inclination [rad]
i = deg2rad(90);

% right ascension of the ascending node [rad]
RAAN = deg2rad(0);

% argument of perigee [rad]
omega = deg2rad(45);

% true anomaly [rad]
nu = deg2rad(180);

oe = [a, e, i, RAAN, omega, nu];
[r_I, v_I] = oe2rv(GM, oe);
X0         = [r_I;v_I];

%% time vector
n_orbits = 1;
dt       = 0.1;   % [sec] 
T        = 2* pi * sqrt(a^3 / GM); % [sec]
t_max    = T * n_orbits;
t_int    = 0:dt:t_max;

%% propagate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13); Nx = 6;
[time, state] = ode113(@(t, x) EOM_2BP(x, t, GM), t_int, X0, options);
r = state(:, 1:3)'; v = state(:, 4:6)';

figure()
plot3(r(1, :), r(2, :), r(3, :));
axis equal; grid on;

%% compute oe time series
oe_time_series = nan(length(time), 6);
for k = 1 :length(time)
    oe_time_series(k, :) = rv2oe(GM, r(:, k), v(:, k));
end

figure();
tt_fig = ["a [m]", "e [-]", "i [rad]",...
    "RAAN [rad]", "\omega [rad]", "\nu [rad]"];
for j = 1:6
    subplot(2, 3, j);
    plot(time, oe_time_series(:, j)); 
    grid on; 
    xlabel('Seconds');
    title(tt_fig(j));
end
sgtitle('Orbital Elements time evolution');

%% compute integrated GG meas. in radial-radial direction
% Integration time [sec]
Tint = 120; 
y_rr_int = nan(1, length(time));
y_rr_ist = nan(1, length(time));
for k = 1:length(time)
    a     = oe_time_series(k, 1);
    e     = oe_time_series(k, 2);
    f0    = oe_time_series(k, 6);
    y_rr_int(k) = int_Grr_2BP(GM, a, e, f0, Tint)./Tint;
    y_rr_ist(k) = meas_Grr_2BP(GM, r(:, k));

end
figure()
subplot(2, 1, 1);
plot(time./T, y_rr_int./1E-9, 'LineWidth', 1.5, 'Color', 'b');
hold on;
plot(time./T, y_rr_ist./1E-9, 'LineWidth', 1.5, 'Color', 'g');
grid on; xlabel('Orbit revolutions'); ylabel('Eotvos');
legend('averaged', 'instantenous');
title('Averaged vs Integrated \Gamma_{rr} over ' + string(Tint) + ' sec.');

subplot(2, 1, 2);
plot(time./T, (y_rr_ist-y_rr_int)./1E-9, 'LineWidth', 1.5, 'Color', 'r');
grid on; xlabel('Orbit revolutions'); ylabel('Eotvos');
title('Averaged and Integrated difference');

%% numerical approximation
Ns_avg = round(Tint/dt);

if abs(Ns_avg*dt - Tint) > 1e-12
    warning(['Tint is not an integer multiple of dt. ' ...
             'Using Ns_avg = %d steps.'], Ns_avg);
end

% RK4 integration step
h_avg = Tint/Ns_avg;

Nt = length(time);

y_rr_rk4   = nan(1, Nt);
y_rr_mean  = nan(1, Nt);
y_rr_ode113 = nan(1, Nt);

% ODE113 tolerances for augmented state:
% [position; velocity; integrated measurement]
options_int = odeset( ...
    'RelTol', 1e-13, ...
    'AbsTol', [1e-6*ones(3,1); ...
               1e-9*ones(3,1); ...
               1e-16]);
% loop
for k = 1:Nt

    % Initial Cartesian state at beginning of averaging window
    Xk = state(k, :)';

    % Integrated measurement using fixed-step RK4
    Irr_rk4 = int_Grr_2BP_RK4(GM, Xk, Tint, Ns_avg);

    y_rr_rk4(k) = Irr_rk4 / Tint;

    % Integrated measurement using ODE113
    Z0 = [Xk; 0];

    [~, Z_ode113] = ode113( ...
        @(tau, Z) EOM_2BP_integrated_measurement(tau, Z, GM), ...
        [0, Tint], Z0, options_int);

    Irr_ode113 = Z_ode113(end, 7);

    y_rr_ode113(k) = Irr_ode113 / Tint;

    % Arithmetic mean of discrete measurements
    % k:k+Ns_avg contains Ns_avg+1 samples, including both endpoints.
    if k + Ns_avg <= Nt
        y_window = y_rr_ist(k:k+Ns_avg);

        y_rr_mean(k) = mean(y_window);
    end
end

% Relative errors
RK4_error = abs(y_rr_int - y_rr_rk4) ./ abs(y_rr_int) * 100;

ODE113_error = abs(y_rr_int - y_rr_ode113) ...
             ./ abs(y_rr_int) * 100;

Mean_error = abs(y_rr_int - y_rr_mean) ...
           ./ abs(y_rr_int) * 100;

% Plot integrated-measurement approximations
figure();

plot(time./T, y_rr_int./1E-9, ...
    'LineWidth', 1.8, 'Color', 'k');
hold on;

plot(time./T, y_rr_rk4./1E-9, ...
    '--', 'LineWidth', 1.4, 'Color', 'g');

plot(time./T, y_rr_ode113./1E-9, ...
    '-.', 'LineWidth', 1.4, 'Color', 'b');

plot(time./T, y_rr_mean./1E-9, ...
    ':', 'LineWidth', 1.6, 'Color', 'r');

grid on;
xlabel('Orbit revolutions');
ylabel('Eotvos');

title(['Integrated \Gamma_{rr} over ', ...
       num2str(Tint), ' sec']);

legend('Analytical', ...
       'RK4', ...
       'ODE113 augmented state', ...
       'Discrete arithmetic mean', ...
       'Location', 'best');

%% Plot relative errors
figure();

semilogy(time./T, RK4_error, ...
    'LineWidth', 1.5, 'Color', 'g');
hold on;

semilogy(time./T, ODE113_error, ...
    'LineWidth', 1.5, 'Color', 'b');

semilogy(time./T, Mean_error, ...
    'LineWidth', 1.5, 'Color', 'r');

grid on;
xlabel('Orbit revolutions');
ylabel('Relative error [%]');
title('Numerical approximation relative error');

legend('RK4', ...
       'ODE113 augmented state', ...
       'Discrete arithmetic mean', ...
       'Location', 'best');

%% helpers
% -------------------------------------------------------------------------
function [dx] = EOM_2BP(x, t, GM)
    % radius vector
    r       = [x(1);x(2);x(3)];
    r_norm = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    
    r_hat = r./r_norm;
    
    % acceleration vector
    a = - GM / (r_norm^2) * r_hat;

    % state time derivative
    dx = [x(4:6);a];
end

function [r_I, v_I] = oe2rv(mu, oe)
    %OE2RV Convert classical orbital elements to Cartesian state in 2BP.
    %
    %   [r_I, v_I] = oe2rv(mu, oe)
    %
    %   Inputs:
    %       mu  - gravitational parameter [km^3/s^2]
    %       oe  - orbital elements:
    %             oe = [a, e, i, RAAN, omega, nu]
    %
    %             a     = semi-major axis [km]
    %             e     = eccentricity [-]
    %             i     = inclination [rad]
    %             RAAN  = right ascension of ascending node [rad]
    %             omega = argument of periapsis [rad]
    %             nu    = true anomaly [rad]
    %
    %   Outputs:
    %       r_I - inertial Cartesian position [km], 3x1
    %       v_I - inertial Cartesian velocity [km/s], 3x1
    
    % Extract orbital elements
    a     = oe(1);
    e     = oe(2);
    inc   = oe(3);
    RAAN  = oe(4);
    omega = oe(5);
    nu    = oe(6);
    
    % Semi-latus rectum
    p = a * (1 - e^2);
    
    % Radius magnitude
    r_mag = p / (1 + e*cos(nu));
    
    % Position and velocity in perifocal frame
    r_PQW = r_mag * [cos(nu);
        sin(nu);
        0];
    
    v_PQW = sqrt(mu/p) * [-sin(nu);
        e + cos(nu);
        0];
    
    % Rotation matrices
    R3_RAAN = [ cos(RAAN), -sin(RAAN), 0;
        sin(RAAN),  cos(RAAN), 0;
        0,          0,         1];
    
    R1_inc = [1, 0,         0;
        0, cos(inc), -sin(inc);
        0, sin(inc),  cos(inc)];
    
    R3_omega = [ cos(omega), -sin(omega), 0;
        sin(omega),  cos(omega), 0;
        0,           0,          1];
    
    % Perifocal to inertial rotation
    Q_PQW_to_I = R3_RAAN * R1_inc * R3_omega;
    
    % Cartesian state in inertial frame
    r_I = Q_PQW_to_I * r_PQW;
    v_I = Q_PQW_to_I * v_PQW;

end

function oe = rv2oe(mu, r_I, v_I)
%RV2OE Convert inertial Cartesian state to classical orbital elements.
%
%   oe = rv2oe(mu, r_I, v_I)
%
%   Inputs:
%       mu  - gravitational parameter [km^3/s^2]
%       r_I - inertial position vector [km], 3x1
%       v_I - inertial velocity vector [km/s], 3x1
%
%   Output:
%       oe  - classical orbital elements:
%             oe = [a, e, i, RAAN, omega, nu]
%
%             a     = semi-major axis [km]
%             e     = eccentricity [-]
%             i     = inclination [rad]
%             RAAN  = right ascension of ascending node [rad]
%             omega = argument of periapsis [rad]
%             nu    = true anomaly [rad]

    tol = 1e-12;

    % Ensure column vectors
    r_I = r_I(:);
    v_I = v_I(:);

    % Magnitudes
    r = norm(r_I);
    v = norm(v_I);

    % Specific angular momentum
    h_vec = cross(r_I, v_I);
    h = norm(h_vec);

    % Inclination
    inc = acos(h_vec(3)/h);

    % Node vector
    k_hat = [0; 0; 1];
    n_vec = cross(k_hat, h_vec);
    n = norm(n_vec);

    % Eccentricity vector
    e_vec = (1/mu) * ((v^2 - mu/r)*r_I - dot(r_I, v_I)*v_I);
    e = norm(e_vec);

    % Specific orbital energy
    energy = v^2/2 - mu/r;

    % Semi-major axis
    if abs(e - 1) > tol
        a = -mu/(2*energy);
    else
        a = Inf; % parabolic case
    end

    % RAAN
    if n > tol
        RAAN = atan2(n_vec(2), n_vec(1));
        RAAN = wrapTo2PiLocal(RAAN);
    else
        RAAN = 0;
    end

    % Argument of periapsis
    if n > tol && e > tol
        omega = atan2( ...
            dot(cross(n_vec, e_vec), h_vec)/h, ...
            dot(n_vec, e_vec));
        omega = wrapTo2PiLocal(omega);
    else
        omega = 0;
    end

    % True anomaly
    if e > tol
        nu = atan2( ...
            dot(cross(e_vec, r_I), h_vec)/h, ...
            dot(e_vec, r_I));
        nu = wrapTo2PiLocal(nu);
    else
        % Circular orbit: use argument of latitude instead
        if n > tol
            nu = atan2( ...
                dot(cross(n_vec, r_I), h_vec)/h, ...
                dot(n_vec, r_I));
            nu = wrapTo2PiLocal(nu);
        else
            % Circular equatorial orbit: use true longitude
            nu = atan2(r_I(2), r_I(1));
            nu = wrapTo2PiLocal(nu);
        end
    end

    % Output orbital elements
    oe = [a, e, inc, RAAN, omega, nu];

end

function ang = wrapTo2PiLocal(ang)
%WRAPTO2PILOCAL Wrap angle to [0, 2*pi).

    ang = mod(ang, 2*pi);

end

function [Irr] = int_Grr_2BP(mu, a, e, f0, Tint)
    %INT_GRR_2BP Analytical time integral of radial-radial gravity gradient.
    %
    %   [Irr, Grr_avg, f1, E1, M1] = int_Grr_2BP(mu, a, e, f0, Tint)
    %
    %   Computes:
    %
    %       Irr = integral Gamma_rr dt
    %
    %   where:
    %
    %       Gamma_rr = 2*mu/r^3
    %
    %   for a Keplerian elliptic orbit.
    %
    %   Inputs:
    %       mu    - gravitational parameter [km^3/s^2]
    %       a     - semi-major axis [km]
    %       e     - eccentricity [-]
    %       f0    - initial true anomaly [rad]
    %       Tint  - integration time [s]
    %
    %   Outputs:
    %       Irr      - integrated rr gravity-gradient measurement [1/s]
    
        if e < 0 || e >= 1
            error('This function assumes an elliptic orbit: 0 <= e < 1.');
        end
    
        if a <= 0
            error('Semi-major axis must be positive for an elliptic orbit.');
        end
    
        if Tint <= 0
            error('Integration time Tint must be positive.');
        end
    
        % Semi-latus rectum
        p = a*(1 - e^2);
    
        % Mean motion
        n = sqrt(mu/a^3);
    
        % Initial eccentric anomaly from true anomaly
        E0 = 2*atan2( ...
            sqrt(1 - e)*sin(f0/2), ...
            sqrt(1 + e)*cos(f0/2));
    
        % Initial mean anomaly
        M0 = E0 - e*sin(E0);
    
        % Final mean anomaly after integration time
        M1 = M0 + n*Tint;
    
        % Solve Kepler equation for final eccentric anomaly
        E1 = solveKeplerElliptic(M1, e);
    
        % Final true anomaly
        f1 = 2*atan2( ...
            sqrt(1 + e)*sin(E1/2), ...
            sqrt(1 - e)*cos(E1/2));
    
        % Make sure delta f is continuous and positive for forward time
        df = f1 - f0;
    
        while df < 0
            df = df + 2*pi;
        end
    
        % If Tint spans more than one orbit, include full revolutions
        Norb = floor((M1 - M0)/(2*pi));
        df_mod = f1 - f0;
    
        while df_mod < 0
            df_mod = df_mod + 2*pi;
        end
    
        df_total = df_mod + 2*pi*Norb;
    
        % Avoid double-counting if df_mod already contains the first revolution
        if df_total > (M1 - M0)/n * n * 2*pi/(2*pi) + 2*pi
            df_total = df;
        end
    
        % Cleaner total true-anomaly increment
        df_total = df_mod + 2*pi*floor((M1 - M0 - df_mod*0)/(2*pi));
    
        % Analytical integral
        Irr = 2*sqrt(mu)/p^(3/2) * ...
              (df_total + e*(sin(f1) - sin(f0)));
end

function E = solveKeplerElliptic(M, e)
    %SOLVEKEPLERELLIPTIC Solve M = E - e sin(E) for elliptic orbits.
    
        % Wrap M for numerical stability
        M_wrapped = mod(M, 2*pi);
    
        % Initial guess
        if e < 0.8
            E = M_wrapped;
        else
            E = pi;
        end
    
        % Newton iterations
        tol = 1e-13;
        maxIter = 50;
    
        for k = 1:maxIter
            F  = E - e*sin(E) - M_wrapped;
            dF = 1 - e*cos(E);
    
            dE = -F/dF;
            E = E + dE;
    
            if abs(dE) < tol
                break;
            end
        end
    
        % Add back full revolutions
        E = E + 2*pi*floor(M/(2*pi));

end

function Grr = meas_Grr_2BP(mu, r_I)
    %MEAS_GRR_2BP Instantaneous radial-radial gravity-gradient measurement.
    %
    %   Grr = meas_Grr_2BP(mu, r_I)
    %
    %   Inputs:
    %       mu  - gravitational parameter [km^3/s^2]
    %       r_I - inertial position vector [km], 3x1
    %
    %   Output:
    %       Grr - radial-radial gravity-gradient component [1/s^2]
    %
    %   Notes:
    %       Assumes 2-body point-mass gravity.
    %       Assumes rr direction is aligned with local radial direction.
    %
    %       Convention used:
    %           Gamma = d(a)/d(r)
    %
    %       Therefore:
    %           Gamma_rr = 2*mu/r^3
    %
    %       If your convention is Gamma = -d(a)/d(r), multiply output by -1.
    
    % Ensure column vector
    r_I = r_I(:);
    
    % Radius magnitude
    r = norm(r_I);
    
    % Radial-radial gravity-gradient measurement
    Grr = 2*mu/r^3;

end

function Irr = int_Grr_2BP_RK4(mu, X0, Tint, Ns)
    %INT_GRR_2BP_RK4 Numerically integrate Gamma_rr over Tint using RK4.
    %
    %   Irr = int_Grr_2BP_RK4(mu, X0, Tint, Ns)
    %
    %   Inputs:
    %       mu   - gravitational parameter [m^3/s^2]
    %       X0   - Cartesian state at beginning of averaging window [r; v]
    %              r in [m], v in [m/s]
    %       Tint - averaging/integration time [s]
    %       Ns   - number of RK4 integration steps inside Tint
    %
    %   Output:
    %       Irr  - integrated radial-radial gravity-gradient measurement [1/s]
    %
    %   The averaged measurement is:
    %
    %       Grr_avg = Irr/Tint
    %
    
    if Tint <= 0
        error('Tint must be positive.');
    end
    
    if Ns < 1
        error('Ns must be at least 1.');
    end
    
    h = Tint/Ns;
    
    % Augmented state:
    % Z = [r; v; I]
    % where I_dot = Gamma_rr
    Z = [X0(:); 0];
    
    for j = 1:Ns
        k1 = EOM_2BP_aug_Grr(mu, Z);
        k2 = EOM_2BP_aug_Grr(mu, Z + 0.5*h*k1);
        k3 = EOM_2BP_aug_Grr(mu, Z + 0.5*h*k2);
        k4 = EOM_2BP_aug_Grr(mu, Z + h*k3);
    
        Z = Z + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end
    
    Irr = Z(7);

end

function dZ = EOM_2BP_aug_Grr(mu, Z)
    %EOM_2BP_AUG_GRR Augmented 2BP dynamics for integrated Gamma_rr.
    %
    %   Z = [r; v; I]
    %
    %   I_dot = Gamma_rr = 2*mu/r^3
    
    r = Z(1:3);
    v = Z(4:6);
    
    r_norm = norm(r);
    
    a = -mu*r/r_norm^3;
    
    Grr = 2*mu/r_norm^3;
    
    dZ = [v; a; Grr];

end

function dZ = EOM_2BP_integrated_measurement(~, Z, GM)
    %% AUGMENTED 2BP DYNAMICS WITH INTEGRATED GRAVITY-GRADIENT MEASUREMENT
    %
    % State:
    %   Z = [r; v; I_rr]
    %
    % where
    %   dI_rr/dt = Gamma_rr
    %
    % and for a point-mass gravity field:
    %   Gamma_rr = 2*GM/r^3
    
    % Extract Cartesian state
    r = Z(1:3);
    v = Z(4:6);
    
    r_norm = norm(r);
    
    % Two-body acceleration
    acceleration = -GM * r / r_norm^3;
    
    % Radial-radial gravity-gradient component
    Gamma_rr = 2 * GM / r_norm^3;
    
    % Augmented state derivative
    dZ = [v;
        acceleration;
        Gamma_rr];
end