clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_navigation/data/')
addpath('../../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/')
addpath('../data/')
set(0,'defaultAxesFontSize',16);

%%            POSITION ERROR EFFECT IN GRAV. ESTIMATION
% Description: Evaluate the position error effect in the gravity estimation
% using a consider covariance formulation
% Author: Sergio Coll
% Date: 01/18/25

% Asteroid parameters.
% % savedData = 0;
% % path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % GM = 5.2;
% % n_max  = 6;
% % normalized = 1;
% % W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
% % W0 = 0;                   % Initial asteroid longitude
% % RA = deg2rad(86.6388);    % Right Ascension     [rad]
% % DEC = deg2rad(-65.1086);  % Declination         [rad]

% Earth parameters
savedData = 1;                % use saved data. 1 = yes / 0 = no
path = "HARMCOEFS_EARTH_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
path = "SIGMACOEFS_EARTH_1.txt";
[sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
GM = 3.986004418E14;
n_max  = 50;
normalized = 1;
W = 2 * pi / (24*3600);     % Rotation ang. vel   [rad/s]
W0 = 0;                     % Initial asteroid longitude
RA = deg2rad(15);           % Right Ascension     [rad]
DEC = 0;                    % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

% Initial conditions
r = Re + 250E3;         % [m] 
r0 = [r;0;0];           % [ACI]
v0 = [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 15*16;
f = 1/10;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% measurement uncertianty
sigmaM = 1E-12;                          % [1/s^2]
noise0 = zeros(9, Nt);

% Integrate trajectory
Earth_case = 0;
if savedData
    load('Nov_L2position.mat');   % ECEF coordinates
    load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF
    rn_ECEF = positions(:, 2:end)';
    t = positions(:, 1);
    Nt = length(t);

    % get rotation matrix
    ACI_ECEF = compute_ECEF2ITRF(outPut);
    rn = rotate2ECI(rn_ECEF, ACI_ECEF, t);
    
    Earth_case = 1;
    t_UTC = convertGPSsc2UTC(t);
else
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    Nx = 6;
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
    [~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 6, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
    rn = state_t(:, 1:3)';
    vn = state_t(:, 4:6)';
    t_UTC = t./86400;   % [days]
end

rn_ECEF = rn.*0;
lat = zeros(1, length(rn(1, :))); lon = lat;
for j = 1:length(rn(1, :))
    if(Earth_case == 0)
        Wt = W0 + W * t(j);
        ACAF_ACI = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    else
        maxPos = 3*j; minPos = maxPos - 2;
        [ACI_ACAF] = ACI_ECEF(minPos:maxPos, :);
        ACAF_ACI = ACI_ACAF'; 
    end
    rn_ECEF(:, j) = ACAF_ACI * rn(:, j);
    lla = ecef2lla(rn_ECEF(:, j)');
    lat(j) = lla(:,1);
    lon(j) = lla(:,2);
end

% Create the 2D plot ground track
figure;
geoplot(lat, lon, '.b') 
geobasemap('landcover')

figure()
plot(t_UTC, (vecnorm(rn) - Re)./1E3, 'LineWidth', 2)
xlabel('Time [days]')
ylabel('[km]')
title('S/C orbital Altitude')

% perturb nominal coefficient
[X_K] = compute_Kaula(n_max, 1E-5);
P0 = diag(X_K(2:end).^2); 
S = X_K;

% attitude error
option = "periodic";    % option: constant / periodic / gaussian
sigmaAt  = 4.84814e-5;    % [rad] 
[At] = compute_attErr(option, t, sigmaAt, T);

% Consider covariance
Pxc = zeros(Ncs-1, 3); Pc = zeros(3, 3);
Pxc_NSM = zeros(Ncs - 1, 6); Pc_NSM = zeros(6, 6);
for j = 1:3
    Pc(j, j) = max(At(1, :))^2;
end
for j = 1:length(Pc_NSM)
    Pc_NSM(j, j) = max(At(1, :))^4;
end

[~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
[~, Mxc_NSM, Mcc_NSM] = get_considerCov_apriori(P0, Pc_NSM, Pxc_NSM);
Ax = inv(P0);  Ax_NSM = inv(P0); sxc = 0;
R0 = diag([sigmaM, sigmaM, sigmaM, sigmaM, sigmaM].^2);
for j = 1:Nt
    fprintf('Loading ... %.2f\n % ', j/Nt * 100);
    % current position
    rn_ACI = rn(:, j);
    vn     = rn.*0;
    
    % ACAF to ACI rotation matrix
    if(Earth_case == 0)
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    else
        maxPos = 3 * j; minPos = maxPos - 2;
        [ACI_ACAF] = ACI_ECEF(minPos:maxPos, :);
        ACAF_ACI = ACI_ACAF'; 
    end

    % computed meas. partials
    [Y, Hx_ACI, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm, eye(3,3));
    hx = [Hx_ACI(1, 2:end); Hx_ACI(2, 2:end); Hx_ACI(3, 2:end);Hx_ACI(5, 2:end);...
        Hx_ACI(8, 2:end)];
    
    % compute attitude partials. % Inertially fixed
    [Hrot_grad] = compute_rotPartials_analy(Y);
    Hrot = Hrot_grad;
    hrot = [Hrot(1, :);Hrot(2,:);Hrot(3,:);Hrot(5, :);Hrot(6, :)];
    
    % select covarinace and estimation parameters
    hc  = hrot; % consider parameters matrix & Attitude NSM
    hcc = zeros(5, 6);

    % LS covariance
    Ax  = Ax  + (hx' * inv(R0) * hx);
    Mxc = Mxc + (hx' * inv(R0) * hc);
    Mcc = Mcc + (hc' * inv(R0) * hc); 

    % NSM covariance
    C = null(hc');
    hx_NSM = C' * hx;   % grav. field partials
    hcc_NSM = C' * hcc;
    r  = C' * R0 * C;
    
    Ax_NSM = Ax_NSM + hx_NSM' * inv(r) * hx_NSM;
    Mxc_NSM = Mxc_NSM + (hx_NSM' * inv(r) * hcc_NSM);
    Mcc_NSM = Mcc_NSM + (hcc_NSM' * inv(r) * hcc_NSM); 

    % error sensitivity
    sxc = sxc + hx' * inv(R0) * (hc * At(:, j));
end

% compute final covariance at epoch time. LS
Px = inv(Ax);
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';
Pxc = Sxc * Pc;
P0_LS = [Pxx, Pxc;Pxc', Pc];

sensitivity = -Px * sxc;

% compute final covariance at epoch time. NSM
Px_NSM = inv(Ax_NSM);
Sxc_NSM = -Px_NSM * Mxc_NSM;
Pxx_NSM = Px_NSM + Sxc_NSM*Pc_NSM*Sxc_NSM';
Pxc_NSM = Sxc_NSM * Pc_NSM;
P0_NSM = [Pxx_NSM, Pxc_NSM;Pxc_NSM', Pc_NSM];

% compute standart deviation (sigma)
s0     = sqrt(diag(P0_LS));
s0_NSM = sqrt(diag(P0_NSM));

% compute RMS value
sigma_RMS_LS  = computeRMS_coeffErr(n_max, Nc, Ns, [1;s0], Cnm.*0, Snm.*0); 
sigma_RMS_NSM = computeRMS_coeffErr(n_max, Nc, Ns, [1;s0_NSM], Cnm.*0, Snm.*0);
sigma_RMS     = computeRMS_coeffErr(n_max, Nc, Ns, S, Cnm.*0, Snm.*0);
RMS = computeRMS_coeffErr(n_max, Nc, Ns, X, Cnm.*0, Snm.*0);
RMS_sensitivity = computeRMS_coeffErr(n_max, Nc, Ns, [1;sensitivity], Cnm.*0, Snm.*0); 

% plot SH RMS value
figure()
semilogy(2:n_max, RMS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','k')
hold all;
semilogy(2:n_max, sigma_RMS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','k', 'LineStyle', '-.')
semilogy(2:n_max, sigma_RMS_LS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','g')
semilogy(2:n_max, sigma_RMS_NSM(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','b')
grid on;
legend('truth', 'apriori', 'LS', 'NSM')
title('RMS value SH coefficients')

% plot sensitivity error ratio
figure()
semilogy(2:n_max, RMS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','k')
hold on;
semilogy(2:n_max, RMS_sensitivity(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','r')
title('Attitude sensitivity')
legend('truth', 'sensitivity')

% plot SH estimation
tt1 = 'Cnm uncertainty';
tt2 = 'Snm uncertainty';
ls  = '-'; mk = 'square'; 
lgn =  {'truth','LS', 'NSM'};
plot_gravField(X, s0, s0_NSM, n_max, tt1, tt2, ls, mk, lgn);


%% FUNCTIONS

function [num_C, num_S, str_C, str_S] = SH_xlabel(n_max)    
    num_C = ones(1, n_max-1) * NaN;
    num_S = num_C;

    str_C = cell(1, n_max - 1);
    str_S = str_C;

    num_C(1) = 1;
    for j = 3:n_max
        num_C(j-1) = j + num_C(j-2);
    end
    
    num_S(1) = 1;
    for j = 3:n_max
        num_S(j-1) = (j-1) + num_S(j-2);
    end

    for j = 2:n_max
        str_C{j - 1} = "C_{" + string(j) + "0}";
        str_S{j - 1} = "S_{" + string(j) + string(j) + "}";
    end 
end

function [] = plot_gravField(X, SH_R, SH_N, n_max, tt1, tt2, ls, mk, lgn)
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    figure()
    subplot(1, 2, 1)
    semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Nc-1, abs(SH_R(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    title(tt1)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    
    subplot(1, 2, 2)
    semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold on;
    semilogy(1:Ns, abs(SH_R(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
    title(tt2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    legend(lgn);
end

function [Hrot] = compute_rotPartials_analy(Y)
    Yxx = Y(1); Yxy = Y(2); Yxz = Y(3);
    Yyy = Y(5); Yyz = Y(6); Yzz = Y(9);

    Hrot = -[0, 2*Yxz, -2*Yxy;...
           -Yxz, Yyz, Yxx - Yyy;...
        Yxy, Yzz - Yxx, -Yyz;...
        -Yxz, Yyz, Yxx - Yyy;...
        -2*Yyz, 0, 2*Yxy;...
        Yyy - Yzz, -Yxy, Yxz;...
        Yxy, Yzz - Yxx, -Yyz;...
        Yyy - Yzz, -Yxy, Yxz;...
        2*Yyz, -2*Yxz, 0];
end

function [At] = compute_attErr(option, t, sigma, T)
    Nt = length(t);
    At = ones(3, Nt) * NaN;
    if(option == "constant")
        for j = 1:Nt
            At(:, j) = sigma;
        end
    elseif(option == "periodic")
        w = 2*pi / T;
        for j = 1:Nt
            At(:, j) = sigma * sin(50*w * t(j));
        end
    elseif(option == "gaussian")
        At = normrnd(0, sigma, [3, Nt]);
    end
end

function [rn_ECI] = rotate2ECI(rn_ECEF, ACI_ECEF, t)
    rn_ECI = ones(3, length(t)) * NaN;

    for j = 1:length(t)
        % rotation matrix
        maxPos = 3*j; minPos = maxPos- 2;
        R = ACI_ECEF(minPos:maxPos, :);

        rn_ECI(:, j) = R * rn_ECEF(:, j);
    end
end

function [t_UTC] = convertGPSsc2UTC(t)
    gps_epoch = datetime(1980,1,6,0,0,0); % GPS epoch
    t_UTC = gps_epoch + seconds(t);
end

function [ACI_ECEF] = compute_ECEF2ITRF(R)
    Nt = length(R(:, 1));
   
    % output matrix
    ACI_ECEF  = ones(3*Nt, 3) * NaN;
    for j = 1:Nt
        b1 = R(j, 2); b2 = R(j, 3); b3 = R(j, 4); b0 = R(j, 5);
       
        % construct matrix
        mat11 = b0^2 + b1^2 - b2^2 - b3^2;
        mat12 = 2*(b1*b2 + b0*b3);
        mat13 = 2*(b1*b3 - b0*b2);
        mat21 = 2*(b1*b2 - b0*b3);
        mat22 = b0^2 - b1^2 + b2^2 - b3^2;
        mat23 = 2*(b2*b3 + b0*b1);
        mat31 = 2*(b1*b3 + b0*b2);
        mat32 = 2*(b2*b3 - b0*b1);
        mat33 = b0^2 - b1^2 - b2^2 + b3^2;

        % ensamble matrix
        maxPos = 3 * j; minPos = maxPos - 2;
        ACI_ECEF(minPos:maxPos, :) = [mat11, mat12, mat13;mat21, mat22, mat23;...
            mat31, mat32, mat33];
    end
end