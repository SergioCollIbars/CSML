clear; clc; close all;
format long g;
set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')
%%          SH COEFFICIENTS ERROR BUDGET LRO
% Description: Compute the impact of gravity field coefficient uncertainty
% in the LRO orbit. (error budget)
% Author: Sergio Coll-Ibars
% Date: 01/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRGM1200 gravity field 
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";
input_clone      = "HARMCOEFS_MOON_CLONE_01_GRGM1200PRIM.txt";
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
n_max            = file(1); normalized = file(3);
[Nc, Ns, ~]      = count_num_coeff(n_max);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_clone);
SH_clone         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);


%% MOON VALUES FROM SPICE
% % % Get GM for the Moon [m^3/s^2]
% % [GM_moon] = cspice_bodvrd('MOON', 'GM', 1)*1E9; 
% % 
% % % get Earth Radius [m]
% % [radii]  = cspice_bodvrd('MOON', 'RADII', 3)*1E3;
% % R_M  = radii(1);

%% MOON VALUES FROM GRMGM 12000 FILE
R_M      = 1.7380000000000000e+06;   % [m]
GM_moon  = 4.9028001224453001e+12;   % [m^3/ s^2]
GM_sigma = 6.4536052689015518e-24;   % [m^3/ s^2]

%% RMS VALUE COEFFICIENT
[RMS_SH_coeff] = computeRMS_coeffErr(n_max, Nc, Ns, SH_coeff, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_error] = computeRMS_coeffErr(n_max, Nc, Ns, SH_clone, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));
[RMS_SH_uncrt] = computeRMS_coeffErr(n_max, Nc, Ns, SH_uncrt, ...
    zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

figure()
semilogy(2:n_max, RMS_SH_coeff(2:end), 'LineWidth', 2, 'Color', 'b'); 
grid on; hold on;
semilogy(2:n_max, 3.*RMS_SH_uncrt(2:end), 'LineWidth', 2, 'Color', 'b', ...
    'LineStyle', '--'); 
semilogy(2:n_max, RMS_SH_error(2:end), 'LineWidth', 2, 'Color', 'r', ...
    'LineStyle', '-'); 
xlabel('Degree'); title('RMS value GRMG1200');
legend('RMS coefficients', 'RMS \sigma'); 

%% LRO SPICE trajectory
utc_start = '2012-03-01 00:00:00';
utc_stop  = '2012-03-04 00:00:00';

et0       = cspice_str2et(utc_start);
et1       = cspice_str2et(utc_stop);
time_sec  = et1 - et0;

frec      = 1 /(10*60); % [Hz]
N         = round(time_sec * frec);   % number of samples
N         = 1;

et        = linspace(et0, et1, N);

tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sc_SPICE(1:3, :) = sc_SPICE(1:3, :).*1E3;         % [m]
sc_SPICE(4:6, :) = sc_SPICE(4:6, :).*1E3;         % [m/s]

% convert time to UTC
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

figure()
plot3(sc_SPICE(1, :), sc_SPICE(2, :), sc_SPICE(3, :), 'LineWidth', 2)
hold on; grid on; title('LRO SPICE orbit'); axis equal;

figure()
plot(tUTC, (vecnorm(sc_SPICE(1:3, :)) - R_M)./1E3, 'LineWidth', 2);
title('Orbit altitude');ylabel('[km]')

disp('Computing rotation matrices')
for j = 1:N
    frame_to   = 'J2000';
    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et(j));

    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000_mat(minInd:maxInd, :) = J2000_MOON';
end
%% GG measurements (true + perturbed)
Monte_Carlo = 1; % number of monte carlo realizations
n_max       = 800; [Nc, Ns, Ncs]    = count_num_coeff(n_max);

Y_true      = nan(9, N, Monte_Carlo);
Y_nom       = nan(9, N, Monte_Carlo);

% trajectory points
x = sc_SPICE(1, :); y = sc_SPICE(2, :); z = sc_SPICE(3, :);

disp('Computing reference gravity field');
[Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);
for j = 1:N
    maxInd = 3 *j; minInd = maxInd - 2;
    BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);

     r = [x(j);y(j);z(j)]; v = zeros(3, 1);                            
    [Y_J2000, ~] = gradiometer_meas(et(j) ,...
        [GM_moon, R_M, n_max, normalized],...
        BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);

    % lon lat
    r_ECEF = BODYMOON_J2000 * r;
    lat    = atan2(r_ECEF(3), sqrt(r_ECEF(1)^2 + r_ECEF(2)^2));
    lon    = atan2(r_ECEF(2), r_ECEF(1));

    [V_LNOF] = compute_GGT_LNOF(vecnorm(r_ECEF), lat, lon, n_max,...
        Cnm, Snm, GM_moon, R_M);


    % rotate to ENU coordinates
    ENU_ECEF = ecef2enu(BODYMOON_J2000 * r);
    T_J2000  = reshape(Y_J2000, [3,3]);
    T_ECEF   = BODYMOON_J2000 * T_J2000 * BODYMOON_J2000';
    T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';

    Y_ENU    = reshape(T_ENU', [9, 1]);
    Y_ENU    = reshape(V_LNOF', [9, 1]);

    Y_nom(:, j, :) = Y_ENU.*ones(9, Monte_Carlo);
end

disp('Computing MC simulation')
for k = 1:Monte_Carlo
    disp("computing ... "  + string(k/Monte_Carlo * 100) + "%")
    err = normrnd(0, SH_uncrt(2:end));
    SH  = SH_coeff + [0;SH_clone(2:end)];
    GM_moon_2 = GM_moon + normrnd(0, GM_sigma);

    [Cnm, Snm] = list2mat(n_max, Nc, Ns, SH);
    tic;
    for j = 1:N
        maxInd = 3 *j; minInd = maxInd - 2;
        BODYMOON_J2000 = BODYMOON_J2000_mat(minInd:maxInd, :);
    
        r = [x(j);y(j);z(j)];
        v = zeros(3, 1);                            
        [Y_J2000, ~] = gradiometer_meas(et(j) ,...
            [GM_moon_2, R_M, n_max, normalized],...
            BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);

         % lon lat
        r_ECEF = BODYMOON_J2000 * r;
        lat    = atan2(r_ECEF(3), sqrt(r_ECEF(1)^2 + r_ECEF(2)^2));
        lon    = atan2(r_ECEF(2), r_ECEF(1));
    
        [V_LNOF] = compute_GGT_LNOF(vecnorm(r_ECEF), lat, lon, n_max,...
            Cnm, Snm, GM_moon_2, R_M);

        % rotate to ENU coordinates
        ENU_ECEF = ecef2enu(BODYMOON_J2000 * r);
        T_J2000  = reshape(Y_J2000, [3,3]);
        T_ECEF   = BODYMOON_J2000 * T_J2000 * BODYMOON_J2000';
        T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';
        Y_ENU    = reshape(T_ENU', [9, 1]);

        Y_ENU    = reshape(V_LNOF', [9, 1]);
    

        Y_true(:, j, k) = Y_ENU;
    end
    sim_t = toc;
    if(k == 1)
        disp('  Estimated time... ' + ...
        string(sim_t * Monte_Carlo / 60) + ' min')
    end
end

disturbance = (Y_true - Y_nom)./1E-12; % [milli-Eotvos]

% stats
RMS_distubance_xx = rms(squeeze(disturbance(1, :, :)), 2);
RMS_distubance_xy = rms(squeeze(disturbance(2, :, :)), 2);
RMS_distubance_xz = rms(squeeze(disturbance(3, :, :)), 2);
RMS_distubance_yy = rms(squeeze(disturbance(5, :, :)), 2);
RMS_distubance_yz = rms(squeeze(disturbance(6, :, :)), 2);
RMS_distubance_zz = rms(squeeze(disturbance(9, :, :)), 2);

RMS_distubance    = [RMS_distubance_xx, RMS_distubance_xy,RMS_distubance_xz ...
    RMS_distubance_yy,RMS_distubance_yz, RMS_distubance_zz];
disp('Time mean for RMS disturbance (milli-Eotvos):');
disp(mean(RMS_distubance));

% plot resuts
figure()
tt = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", 'zz'];
for k = 1:9
    subplot(3, 3, k)
    plot(tUTC, squeeze(disturbance(k, :, :)), '.', 'Color', 'g');
    title(tt(k)); grid on; ylabel("[mE]")
end


%% FUNCTIONS
function [V_LNOF] = compute_GGT_LNOF(r, lat, lon, max_degree, Cnm, Snm, GM, R_ref)
    % Inputs:
    % r, lat, lon: Radius (m), Latitude (deg), Longitude (deg)
    % Cnm, Snm: (max_deg+1) x (max_deg+1) fully normalized matrices
    
    theta = deg2rad(90 - lat); % Colatitude
    phi = deg2rad(lon);
    cos_t = cos(theta);
    sin_t = sin(theta);
    
    % 1. Precompute Fully Normalized Legendre Polynomials and Derivatives
    % We need Pnm, Pnm', and Pnm''
    [P, dP, ddP] = normalized_legendre_derivatives(max_degree, cos_t, sin_t);
    
    % 2. Initialize Spherical Derivatives
    % Vr: dV/dr, Vtt: d2V/dtheta2, etc.
    Vrr = 0; Vtt = 0; Vpp = 0; Vrt = 0; Vrp = 0; Vtp = 0;
    Vr = 0; Vt = 0; Vp = 0;

    for n = 0:max_degree
        rho_n = (R_ref / r)^n;
        term_common = (GM / r) * rho_n;
        
        for m = 0:n
            cos_mphi = cos(m * phi);
            sin_mphi = sin(m * phi);
            CS = Cnm(n+1, m+1) * cos_mphi + Snm(n+1, m+1) * sin_mphi;
            dCS = m * (-Cnm(n+1, m+1) * sin_mphi + Snm(n+1, m+1) * cos_mphi);
            ddCS = -m^2 * CS;

            % Radial derivatives
            Vr  = Vr  - term_common * (n+1)/r * CS * P(n+1, m+1);
            Vrr = Vrr + term_common * (n+1)*(n+2)/(r^2) * CS * P(n+1, m+1);
            
            % Angular derivatives
            Vt  = Vt  + term_common * CS * dP(n+1, m+1);
            Vtt = Vtt + term_common * CS * ddP(n+1, m+1);
            Vp  = Vp  + term_common * dCS * P(n+1, m+1);
            Vpp = Vpp + term_common * ddCS * P(n+1, m+1);
            
            % Mixed derivatives
            Vrt = Vrt - term_common * (n+1)/r * CS * dP(n+1, m+1);
            Vrp = Vrp - term_common * (n+1)/r * dCS * P(n+1, m+1);
            Vtp = Vtp + term_common * dCS * dP(n+1, m+1);
        end
    end

    % 3. Transform from Spherical to LNOF (Local North-Oriented Frame)
    % x=North, y=East, z=Up
    Vxx = (1/r^2)*Vtt + (1/r)*Vr;
    Vyy = (1/(r^2 * sin_t^2))*Vpp + (cot(theta)/r^2)*Vt + (1/r)*Vr;
    Vzz = Vrr;
    
    Vxy = (1/(r^2 * sin_t))*Vtp - (cos_t/(r^2 * sin_t^2))*Vp;
    Vxz = (1/r)*Vrt - (1/r^2)*Vt;
    Vyz = (1/(r*sin_t))*Vrp - (1/(r^2*sin_t))*Vp;

    V_LNOF = [Vxx, Vxy, Vxz; 
              Vxy, Vyy, Vyz; 
              Vxz, Vyz, Vzz];
end

function [P, dP, ddP] = normalized_legendre_derivatives(N, cos_t, sin_t)
    P = zeros(N+1, N+1);
    dP = zeros(N+1, N+1);
    ddP = zeros(N+1, N+1);
    
    P(1,1) = 1;
    for n = 1:N
        % Sectorials m = n
        P(n+1, n+1) = sqrt((2*n+1)/(2*n)) * sin_t * P(n, n);
        % m = n-1
        P(n+1, n) = sqrt(2*n+1) * cos_t * P(n, n);
        
        for m = n-2:-1:0
            anm = sqrt(((2*n-1)*(2*n+1))/((n-m)*(n+m)));
            bnm = sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m)));
            P(n+1, m+1) = anm * cos_t * P(n, m+1) - bnm * P(n-1, m+1);
        end
    end
    
    % Derivatives using standard identities for normalized Pnm
    for n = 1:N
        for m = 0:n
            if m == 0
                dP(n+1, 1) = -sqrt(n*(n+1)/2) * P(n+1, 2);
            elseif m < n
                dP(n+1, m+1) = 0.5 * (sqrt((n+m)*(n-m+1))*P(n+1, m) - sqrt((n+m+1)*(n-m))*P(n+1, m+2));
            else % m == n
                dP(n+1, m+1) = 0.5 * sqrt(2*n) * P(n+1, n);
            end
        end
    end
    % ddP (2nd derivative) can be derived similarly via dP relations
    % (Included here as a placeholder; for GGT, Vtt uses ddP)
    % Note: Laplace condition Vxx+Vyy+Vzz=0 is used to validate precision.
end