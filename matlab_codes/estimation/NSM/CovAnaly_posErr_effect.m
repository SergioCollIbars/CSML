clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
addpath('../../../QGG_navigation/data/')
addpath('data/')
set(0,'defaultAxesFontSize',16);

%%            POSITION ERROR EFFECT IN GRAV. ESTIMATION
% Description: Evaluate the position error effect in the gravity estimation
% using a consider covariance formulation
% Author: Sergio Coll
% Date: 01/18/25


% Asteroid parameters.
savedData = 0;
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % % Earth parameters
% % savedData = 0;                % use saved data. 1 = yes / 0 = no
% % path = "HARMCOEFS_EARTH_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % path = "SIGMACOEFS_EARTH_1.txt";
% % [sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
% % GM = 3.986004418E14;
% % n_max  = 120;
% % normalized = 1;
% % W = 2 * pi / (24*3600);     % Rotation ang. vel   [rad/s]
% % W0 = 0;                     % Initial asteroid longitude
% % RA = -pi/2;                 % Right Ascension     [rad]
% % DEC = pi/2;                 % Declination         [rad]

% % % Eros parameters
% % savedData = 0;
% % path = "HARMCOEFS_EROS_CD_1.txt";
% % name = "EROS";
% % [Cnm, Snm, Re] = readCoeff(path);
% % n_max  = 10;
% % normalized = 1;
% % GM =  459604.431484721;          % Point mass value    [m^3/s^2]
% % W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
% % W0 = 0;                          % Initial asteroid longitude
% % RA = deg2rad(11.363);            % Right Ascension     [rad]
% % DEC = deg2rad(17.232);           % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

% Initial conditions
r      = 2E3;         % [m]
% % r = 24E3;               % [m]
% % r = Re + 250E3;         % [m] 
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/30;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% measurement uncertianty
sigma = 3E-12;                          % [1/s^2]
noise0 = zeros(9, Nt);

% Integrate trajectory
if savedData
    rn = load('position_n85_f10_rev100.mat').rn;
    vn = load('velocity_n85_f10_rev100.mat').vn;
else
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
% %     Nx = Ncs + 5;
    Nx = 6;
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
    [~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
    rn = state_t(:, 1:3)';
    vn = state_t(:, 4:6)';
end

rn_ECEF = rn.*0;
lat = zeros(1, length(rn(1, :))); lon = lat;
for j = 1:length(rn(1, :))
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        rn_ECEF(:, j) = ACAF_ACI * rn(:, j);

        racaf = rn_ECEF(:, j) / vecnorm(rn_ECEF(:, j));
        lat(j) = atan2(racaf(3), sqrt(racaf(1)^2 + racaf(2)^2));
        lon(j) = atan2(racaf(2), racaf(1));
end

% Create the 2D plot
figure;
scatter(rad2deg(lon), rad2deg(lat), 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('2D Plot of Latitude and Longitude');
grid on;
hold on;

figure()
plot(t./86400, (vecnorm(state_t(:, 1:3)') - Re)./1E3, 'LineWidth', 2)
xlabel('Time [days]')
ylabel('[km]')
title('S/C orbital Altitude')

% perturb nominal coefficient
sigma_n = 1 * ones(1, n_max);
[~, Pp] = perturb_coeff(sigma_n, n_max, X);
P0 = Pp(2:end, 2:end); 
S = sqrt(diag(Pp));

% % [S] = mat2list(sigma_Cnm, sigma_Snm, Nc, Ns);
% % S = 1.*S;
% % P0 = diag((S(2:end)).^2);

% Consider covariance
Ar = 1.5E-4;   % [m]
Pxc = zeros(Ncs-1, 3); Pc = zeros(3, 3);
Pxc_NSM = zeros(Ncs - 1, 6); Pc_NSM = zeros(6, 6);
for j = 1:3
    Pc(j, j) = (Ar(1))^2;
end
for j = 1:length(Pc_NSM)
    Pc_NSM(j, j) = Ar(1)^4;
end
c = Ar.*1; % apriori values for the Consider Parameters; 
c_NSM = ones(6, 1).*Ar(1)^2.*1;

[~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
[~, Mxc_NSM, Mcc_NSM] = get_considerCov_apriori(P0, Pc_NSM, Pxc_NSM);
Ax = inv(P0);  Ax_NSM = inv(P0);
R0 = diag([sigma, sigma, sigma, sigma, sigma].^2);
for j = 1:Nt
    fprintf('Loading ... %.2f\n % ', j/Nt * 100);
    % current position
    rn_ACI = rn(:, j);
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % computed meas. partials
    [~, Hx_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm, eye(3,3));
    hx = [Hx_ACI(1, 2:end); Hx_ACI(2, 2:end); Hx_ACI(3, 2:end);Hx_ACI(5, 2:end);...
        Hx_ACI(6, 2:end)];

    % compute Consider Params partials for LS and NSM
    [Hc] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, eye(3,3));
    hpos = [Hc(1, :);Hc(2,:);Hc(3,:);Hc(5, :);Hc(6, :)];
    hap = compute_posPartials_2ndOrder(GM, rn_ACI(1), rn_ACI(2), rn_ACI(3));
    
    % compute attitude partials. % Inertially fixed
    [Hrot_grad] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, ACAF_ACI);
    Hrot = Hrot_grad;
    hrot = [Hrot(1, :);Hrot(2,:);Hrot(3,:);Hrot(5, :);Hrot(6, :)];
    
    % select covarinace and estimation parameters
% %     hc = hpos; % consider parameters matrix. Position NSM
    hc = hrot; % consider parameters matrix. Attitude NSM
    hap = hap.*0;

    % LS covariance
    Ax  = Ax  + (hx' * inv(R0) * hx);
    Mxc = Mxc + (hx' * inv(R0) * hc);
    Mcc = Mcc + (hc' * inv(R0) * hc); 

    % NSM covariance
    C = null(hc');
    hx_NSM = C' * hx;
    hap_NSM = C' * hap;
    r  = C' * R0 * C;
    
    Ax_NSM = Ax_NSM + hx_NSM' * inv(r) * hx_NSM;
    Mxc_NSM = Mxc_NSM + (hx_NSM' * inv(r) * hap_NSM);
    Mcc_NSM = Mcc_NSM + (hap_NSM' * inv(r) * hap_NSM); 
end

% compute final covariance at epoch time. LS
Px = inv(Ax);
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';
Pxc = Sxc * Pc;
P0_LS = [Pxx, Pxc;Pxc', Pc];

% compute final covariance at epoch time. NSM
Px_NSM = inv(Ax_NSM);
Sxc_NSM = -Px_NSM * Mxc_NSM;
Pxx_NSM = Px_NSM + Sxc_NSM*Pc_NSM*Sxc_NSM';
Pxc_NSM = Sxc_NSM * Pc_NSM;
P0_NSM = [Pxx_NSM, Pxc_NSM;Pxc_NSM', Pc_NSM];

% compute standart deviation (sigma)
s0     = sqrt(diag(P0_LS));
s0_NSM = sqrt(diag(P0_NSM));

% compute analytical sigma pos value
C = diag(Px_NSM - Px);
A = diag(Sxc_NSM * Sxc_NSM');
B = diag(Sxc * Sxc');
sigmaTh1 = A.*0; sigmaTh2 = sigmaTh1; sigmaTh3 = sigmaTh1;
for j = 1:length(A)
    c = C(j);
    b = B(j);
    a = A(j);
    y = (- b + sqrt(b^2 - 4 * a * c)) / (-2*a);
    if(y < 0), y = (- b - sqrt(b^2 - 4 * a * c)) / (-2*a); end
    
    % option including 'a' and 'b'
    sigmaTh3(j) = sqrt(y);
   
    % approximation including only 'b'
% %     sigmaTh1(j) = (1/sigma)*sqrt(c/b) * 1E-9;                    % [m / E]
    sigmaTh1(j) = (1/sigma)*sqrt(c/b) * 1E-9 * 3600 * 180 / pi;  % [arcSec / E]
    sigmaTh2(j) = sqrt(c/b);                                     % [m]

end
[C_Perr, S_Perr] = list2mat(n_max, Nc, Ns, [0; sigmaTh1]);
C_Perr(1, :) = C_Perr(1, :).*NaN; C_Perr(2, :) = C_Perr(2, :).*NaN;
S_Perr(:, 1) = S_Perr(:, 1).*NaN; S_Perr(1:2, :) = S_Perr(1:2, :).*NaN;

% Plot the values (units of m/Eotvos)
[zonal, sectoral] = get_ZonalSectoral(n_max, C_Perr, S_Perr);

figure()
semilogy(2:n_max, zonal, 'Marker','square', 'LineStyle','--', "MarkerFaceColor", 'auto', ...
    'MarkerEdgeColor','auto', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'r')
grid on;
xticks(2:n_max);
xticklabels(string(2:n_max));
xlabel('Degree')
ylabel('[m / Eotvos]')
title('Zonal coefficients')

figure()
subplot(1, 2, 1)
semilogy(2:n_max, sectoral(1, :), 'Marker','square', 'LineStyle','--', "MarkerFaceColor", 'auto', ...
    'MarkerEdgeColor','auto', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'r')
grid on;
xticks(2:n_max);
xticklabels(string(2:n_max));
xlabel('Degree')
ylabel('[m / Eotvos]')
title('C_{nm} sectoral coefficients')

subplot(1, 2, 2)
semilogy(2:n_max, sectoral(2, :), 'Marker','square', 'LineStyle','--', "MarkerFaceColor", 'auto', ...
    'MarkerEdgeColor','auto', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'r')
grid on;
xticks(2:n_max);
xticklabels(string(2:n_max));
xlabel('Degree')
ylabel('[m / Eotvos]')
title('S_{nm} sectoral coefficients')

% Ensure 0 values appear white (units of m)
[C_Perr, S_Perr] = list2mat(n_max, Nc, Ns, [0; sigmaTh2]);
C_Perr(1, :) = C_Perr(1, :).*NaN; C_Perr(2, :) = C_Perr(2, :).*NaN;
S_Perr(:, 1) = S_Perr(:, 1).*NaN; S_Perr(1:2, :) = S_Perr(1:2, :).*NaN;

C_Perr(C_Perr == 0) = NaN;   S_Perr(S_Perr == 0) = NaN; 

C_Perr(C_Perr < 1e-3) = 1e-3; % Set anything below 1 mm to 1 mm
C_Perr(C_Perr > 100) = 100;   % Set anything above 10 m to 10 m

S_Perr(S_Perr < 1e-3) = 1e-3; % Set anything below 1 mm to 1 mm
S_Perr(S_Perr > 100) = 100;   % Set anything above 10 m to 10 m

figure;
subplot(1, 2, 1)
imagesc(C_Perr);
set(gca, 'ColorScale', 'log'); % Logarithmic scale
colormap("turbo"); % Use a smooth colormap
colorbar;
xticks(linspace(0, n_max, 7)); % Set 6 x-ticks
yticks(linspace(0, n_max, 7))
c = colorbar;
c.Ticks = [1e-3, 1e-2, 1e-1, 1, 10, 100]; % 1mm, 1cm, 10cm, 1m, 10m
c.TickLabels = {'1 mm', '1 cm', '10 cm', '1 m', '10 m', '100m'};
title('C_{nm} coefficients')
xlabel('Degree')
ylabel('Order')

subplot(1, 2, 2)
imagesc(S_Perr);
set(gca, 'ColorScale', 'log'); % Logarithmic scale
colormap("turbo"); % Use a smooth colormap
colorbar;
xticks(linspace(0, n_max, 7)); % Set 6 x-ticks
yticks(linspace(0, n_max, 7))
c = colorbar;
c.Ticks = [1e-3, 1e-2, 1e-1, 1, 10, 100]; % 1mm, 1cm, 10cm, 1m, 10m
c.TickLabels = {'1 mm', '1 cm', '10 cm', '1 m', '10 m', '100 m'};
title('S_{nm} coefficients')
xlabel('Degree')
ylabel('Oder')

figure;
J = ones(n_max+1, (n_max+1)*2) * NaN; B = fliplr(S_Perr);
J(:, 1:n_max+1) = B;
J(:, n_max+2:end) = C_Perr;
im = imagesc(J);
set(gca, 'ColorScale', 'log');
colormap("turbo");
colorbar
yticks(linspace(1, n_max + 1, 7))
yticklabels(compose('%i', linspace(0, 120, 7)));
xticks(linspace(1, 2*n_max + 2, 7))
xticklabels(compose('%i', linspace(-120, 120, 7)));
c = colorbar;
c.Ticks = [1e-3, 1e-2, 1e-1, 1, 10, 100]; % 1mm, 1cm, 10cm, 1m, 10m
c.TickLabels = {'1 mm', '1 cm', '10 cm', '1 m', '10 m', '100m'};
xlabel('Degree');
ylabel('Order');

% compute RMS value
sigma_RMS_LS  = computeRMS_coeffErr(n_max, Nc, Ns, [1;s0], Cnm.*0, Snm.*0); 
sigma_RMS_NSM = computeRMS_coeffErr(n_max, Nc, Ns, [1;s0_NSM], Cnm.*0, Snm.*0);
sigma_RMS     = computeRMS_coeffErr(n_max, Nc, Ns, S, Cnm.*0, Snm.*0);
RMS = computeRMS_coeffErr(n_max, Nc, Ns, X, Cnm.*0, Snm.*0);

figure()
semilogy(2:n_max, RMS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','k')
hold all;
 semilogy(2:n_max, sigma_RMS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','k', 'LineStyle', '--')
semilogy(2:n_max, sigma_RMS_LS(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','g')
semilogy(2:n_max, sigma_RMS_NSM(2:end), 'LineWidth', 2, 'Marker', 'square', 'Color','b')
grid on;
legend('truth', 'apriori', 'LS', 'NSM')
title('RMS value SH coefficients')

% plot SH estimation
tt1 = 'Cnm uncertainty';
tt2 = 'Snm uncertainty';
ls  = '-'; mk = 'square'; 
lgn =  {'truth','LS', 'NSM'};
plot_gravField(X, s0, s0_NSM, n_max, tt1, tt2, ls, mk, lgn);


%% FUNCTIONS
function [Xp, P0] = perturb_coeff(sigma_n, n_max, X)
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 

    % perturbed values
    Xp = X;
    P0 = eye(Ncs);

    m  = 0;
    n = 2;
    for j =2:Nc
        Xp(j) = Xp(j) + normrnd(0, sigma_n(n-1));
        P0(j, j) = sigma_n(n-1)^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        Xp(j) = Xp(j) + normrnd(0, sigma_n(n-1));
        P0(j, j) = sigma_n(n-1)^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

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

function [sigmaN] = computeRMS_analytical(GM, Re, n, w, L, r)
    % define RMS sigma
    sigmaN = ones(1, n) * NaN;
    sigmaN(1) = r^3/(GM) * w/sqrt(L) * sqrt(1/5);

    for j = 2:n
        P = 3*j^4 + 14*j^3 + 47/2*j^2 + 35/2*j + 5;
        sigmaN(j) = r^(3+j)/(GM * Re^j) * w / sqrt(L) * sqrt(1 / P);
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

function [zonal, sectoral] = get_ZonalSectoral(n_max, C, S)
    % output matrices
    zonal  = ones(1, n_max - 1) * NaN;      % from n = 2 to n_max
    sectoral = ones(2, n_max - 1) * NaN;    % from n = 2 to n_max
    
    % fill up zonal coeffcients
    zonal = C(3:end, 1);
    
    % fill up sectoral coefficients
    for j = 2:n_max
        val = C(j+1, j+1);
        if(isreal(val)), sectoral(1, j-1) = val; else, sectoral(1, j-1) = NaN; end

        val = S(j + 1, j + 1);
        if(isreal(val)), sectoral(2, j-1) = val; else, sectoral(2, j-1) = NaN; end
    end
end
