clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');
%% SIGNAL TO NOISE RATIO OF QUANTUM GRADIOMETERS IN EARTH ORBIT
% Author: Sergio Coll-Ibars
% Date: 07/28/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Planet constants

GM_E = 3.986E14; % [m^3/s^2]
R_E  = 6371E3;  % [m]
path_Earth = "/Users/sergiocollibars/Desktop/" + ...
    "CSML/QGG_navigation/LRO_navigation/data/" + ...
    "HARMCOEFS_EARTH_1.txt";
[Cnm_E, Snm_E, ~] = readCoeff(path_Earth);

%% Instrument parameters
Ts = 10;

%% Orbit parameters

% orbit altitude
H   = 300E3;    % [m]

%% Compute SNR
[sigma, ~] = noise_model_constant(Ts);

[y_avrg, time, r, v] = compute_SNR(GM_E, R_E, Cnm_E, Snm_E, H, Ts);
SNR                  = abs(y_avrg)./sigma;

% plot orbit
plot3(r(1, :), r(2, :), r(3, :), 'Color', 'b');
axis equal; grid on;
xlabel('X[Km]'); ylabel('Y[Km]'); zlabel('Z[Km]');

% plot SNR
figure();
for k = 1:6
    subplot(2, 3, k);
    semilogy(time(2:end), SNR(k, :), 'LineWidth', 2, 'Color','g');
    grid on; xlabel('sec')
end

%% Sampling times SCAN per SH DEGREE
Ts    = 0:10:60;
Ts(1) = 1;

componentNames = {'XX','YY','ZZ','XY','XZ','YZ'};

PSD  = cell(length(Ts),6);
Freq = cell(length(Ts),1);

for j = 1:length(Ts)
    disp(Ts(j))
    [sigma, PSD_noise] = noise_model_constant(Ts(j));

    [y_avrg,time,r,v] = compute_SNR( ...
        GM_E,R_E,Cnm_E,Snm_E,H, Ts(j));

    dt = median(diff(time));      % seconds
    fs = 1/dt;

    % Mean orbital speed
    vMean = mean(vecnorm(v,2,1));

    % Frequency vector only needs to be stored once
    for k = 1:6

        x = detrend(y_avrg(k,:),0);

        [Pxx,f] = pwelch(x,[],[],[],fs);

        PSD{j,k} = Pxx./(PSD_noise);
    end

    Freq{j} = f;

    % Approximate SH degree corresponding to this frequency vector
    Degree{j} = pi*R_E*f/vMean;

end

for k = 1:6

    figure;
    ax = axes;
    hold(ax, 'on');

    for j = 1:length(Ts)

        x = Freq{j};       % Or Degree{j}, if already converted to SH degree
        y = PSD{j,k};

        % Ensure column vectors
        x = x(:);
        y = y(:);

        % Remove DC, zero/negative values, NaNs, and Inf values
        valid = x > 0 & y > 0 & isfinite(x) & isfinite(y);

        loglog(ax, x(valid), y(valid), ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('T_s = %d s', Ts(j)));

    end

    ax.XScale = 'log';
    ax.YScale = 'log';

    xlabel(ax, 'Frequency [Hz]')
    ylabel(ax, 'PSD')
    title(ax, componentNames{k})

    grid(ax, 'on')
    legend(ax, 'show', 'Location', 'best')

end


%% AUXILIARY FUNCTIONS
function [y_avrg, time, r, v] = compute_SNR(GM_E, R_E, Cnm_E, Snm_E, H, Ts)
    % orbit initial parameters
    r0  = [R_E+H;0;0];                      % [m]
    v0  = [0;0;sqrt(GM_E / vecnorm(r0))];   % [m/s]
    X0  = [r0;v0;zeros(6,1)];
    
    % time vector
    T   = 2*pi*((R_E+H)^3/GM_E)^(0.5);  % [s]
    % % T   = 3 * 86400;                   % [s]
    N   = floor(T/Ts);
    
    time = (0:N)*Ts;
    
    %% Integrated state + measurement
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~, state] = ode113(@(t, x) EOM_Earth(x, t, GM_E, R_E, Cnm_E, Snm_E),...
        time, X0, options);
    r      = state(:, 1:3)'; v = state(:, 4:6)'; y_acum = state(:, 7:end)';
    y_avrg = diff(y_acum, 1, 2)./Ts;
end

function [dx] = EOM_Earth(x, t, GM, R, Cnm, Snm)
    % state: pos  3x1
    %        vel  3x1
    %        meas 6x1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract state vector
    r_J2000      = x(1:3);
    
    % Earth orientation
    J2000_EARTH = cspice_pxform('IAU_EARTH', 'J2000', t);

    % LNOF orientation
    k_J   = J2000_EARTH * [0; 0; 1];
    C_L_J = J2000_to_LNOF(r_J2000, k_J);
    
    % Earth central gravity
    [~, a_IAU, T_IAU] = potentialGradient_nm_mex( ...
        Cnm, Snm, 120, ...
        J2000_EARTH' * r_J2000, R, GM, 1);
    
    % acceleration to J2000
    a_J2000 = J2000_EARTH * a_IAU;
    
    % translate tensor to J2000
    T_J2000   =  J2000_EARTH * (T_IAU) * J2000_EARTH';

    % rotate to instrument frame
    T_B       = C_L_J * T_J2000 * C_L_J';
    y_vec     = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);...
                 T_B(2,3);T_B(3,3)];
    
    % state time derivative
    dx = [x(4:6);a_J2000;y_vec];
end

function C_L_J = J2000_to_LNOF(r_J, k_J)
    % J2000_to_LNOF
    %
    % Computes the transformation from J2000 to the Local North-Oriented
    % Frame (LNOF).
    %
    % LNOF convention:
    %   x-axis: North
    %   y-axis: West
    %   z-axis: Up
    %
    % Inputs:
    %   r_J : spacecraft position in J2000, 3x1
    %   k_J : central-body north-pole unit vector in J2000, 3x1
    %
    % Output:
    %   C_L_J : DCM mapping J2000 components into LNOF components
    %
    %       v_L = C_L_J * v_J
    
    r_J = r_J(:);
    k_J = k_J(:);
    
    % Normalize central-body pole
    k_J = k_J / norm(k_J);
    
    % Local Up
    u_J = r_J / norm(r_J);
    
    % Local North: pole projected onto tangent plane
    n_J = k_J - dot(k_J, u_J) * u_J;
    
    if norm(n_J) < 1e-12
        error('LNOF is singular near the geographic poles.');
    end
    
    n_J = n_J / norm(n_J);
    
    % Local West
    w_J = cross(u_J, n_J);
    w_J = w_J / norm(w_J);
    
    % Recompute North to improve numerical orthogonality
    n_J = cross(w_J, u_J);
    n_J = n_J / norm(n_J);
    
    % J2000 to LNOF
    C_L_J = [n_J.';
        w_J.';
        u_J.'];
end

function [sigma, PSD] = noise_model_constant(Ts)
    ASD  = 1E-12; 
    sigma = ASD * sqrt(1/Ts); % [s^-2] 

    PSD   = ASD^2 * (1/Ts)^6;
end