clear;
clc;
close all;

% ---- Body constants (Moon defaults; change as needed) ----
J2  = -0.00020322865242891;     % [-]
R   = 1737.4;                   % [km] 
mu  = 4.902800118E3;            % [km^3/s^2] 

% load gravity coefficients
path = "HARMCOEFS_MOON_1_v2.txt";
[Cmat, Smat, ~] = readCoeff(path); % grav. field primary

% gravity field degree & order
n = 2;
m = 0;

% ---- Choose ONE set of orbital elements (radians) ----
a     = R + 50;                  % [km]
e     = 0;                       % set 0 for circular test
i     = deg2rad(90);             % inclination [rad]
Omega = deg2rad(0);              % RAAN [rad]
omega = deg2rad(0);              % arg. perigee [rad]
M     = deg2rad(0);              % mean anomaly [rad]

% ---- (Optional) body sidereal angle to go ECI -> body-fixed (for lambda)
theta = deg2rad(0);               % set to 0 if you like

% ---- Solve Kepler, get f, r, u ----
E  = keplerE(M, e);
f  = 2*atan2( sqrt(1+e).*sin(E/2), sqrt(1-e).*cos(E/2) );
r  = a*(1 - e*cos(E));
u  = omega + f;                   % argument of latitude

% Build ECI position from elements, rotate to body-fixed to get lambda
[rECI, vECI] = coe2rv(a,e,i,Omega,omega,M,mu);

T = 2*pi / sqrt(mu / a^3);        % [seconds]
fs = 1/1;                         % [Hz]
tmax = 500*T;
spacing = 1/fs;

N = round(tmax/spacing + 1);
if(mod(N,2) ~= 0)
    N = N + 1;
end
time = linspace(0, tmax, N);     % [seconds]

% integrate trajectory (point mass gravity field)
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, state] = ode113(@(t, x) integrator_EOM(x, t, mu), time, [rECI;vECI], ...
    options);

figure()
plot3(state(:, 1), state(:, 2), state(:, 3), 'LineWidth',2);
axis equal;
grid on;

figure()
plot(t, vecnorm(state(:, 1:3)'), 'LineWidth', 2);

% compute measurements in rr (J2 component only)
Trr = ones(1, N) * NaN;
for k = 1:N
    rmag = vecnorm(state(k, 1:3)');
    phi  = atan2(state(k, 3), sqrt(state(k, 1)^2 + state(k, 2)^2)); 
    lambda = atan2(state(k, 2), state(k, 1));
    Anm = (n+1)*(n+2) * GM/rmag * (R/rmag)^n;
    Pnm_array = legendre(n,sin(phi));
    Pnm = Pnm_array(m+1);
    Trr(k) = 12*(mu*J2*R^2)/(2*rmag^5) * ( 3*(sin(phi)^2) - 1 );
% %     Trr(k) = (mu*J2*R^2)/(2*rmag^3) * (3*sin(phi)^2 - 1);
end

figure()
plot(time, Trr./1E-9, 'LineWidth', 2);
xlabel('Time [sec]');
ylabel('[Eotvos]');

% compute FFT for the signal
Fs = 1 / (time(2) - time(1));
compute_FFT(Fs, N, Trr./1E-9)

% % [f, A] = plot_fft_ts(fs, Trr.*1E6, struct('window','hann','title','FFT of x'));
% % set(gca, 'XScale', 'log') ;

%% FUNCTIONS
function [dx] = integrator_EOM(x, t, GM)
    r = x(1:3);
    v = x(4:6);

    % compute the gravity acceleration
    a = -GM* r/(vecnorm(r)^3);

    dx = [v; a];
end

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

function [] = compute_FFT(Fs, L, signal)
    % Description: Fs = sampling frequency
    %              L  = number of samples 
    %              signal = time series of the signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x = signal(:) - mean(signal);           % remove DC
    L = length(x);
    w = hann(L,'periodic');                 % window to control leakage
    CG = mean(w);                           % coherent gain (Hann => 0.5)
    xw = x .* w;

    N = L;                                  % use L for amplitude scaling
    X = fft(xw, N);

    f  = Fs*(0:(N/2))/N;
    P2 = abs(X) / (L*CG);                   % correct for window loss
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);            % single-sided

    figure
    plot(f, P1, 'LineWidth', 2);
    xlabel('f (Hz)'); ylabel('Amplitude');
    title('Single-Sided Amplitude Spectrum (window-corrected)');
    set(gca,'XScale','log'); set(gca,'YScale','log'); grid on;
end
