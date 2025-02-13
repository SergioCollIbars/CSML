clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Compare Rummel's formulation vs the null space approach.
% Test in small radius orbits using apriori information
% Author: Sergio Coll
% Date: 09/28/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

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


% % path = "HARMCOEFS_EARTH_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % path = "SIGMACOEFS_EARTH_1.txt";
% % [sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
% % GM = 3.986004418E14;
% % n_max  = 10;
% % normalized = 1;
% % W = 2 * pi / (24*3600);     % Rotation ang. vel   [rad/s]
% % W0 = 0;                     % Initial asteroid longitude
% % RA = -pi/2;                 % Right Ascension     [rad]
% % DEC = pi/2;                 % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.3E3;
% % r      = Re + 300E3; 
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 2e-2*[1;1;1].*1;            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% noise values from GOCE mission
noise0 = zeros(9, Nt);
% % sigma1  = 0.01 * 1E-9 * sqrt(f); % Vxx, Vyy
% % sigma2  = 0.6  * 1E-9 * sqrt(f); % Vyz, Vyx
% % sigma3  = 0.02 * 1E-9 * sqrt(f); % Vxz, Vzz

sigma1 = 1E-15;
sigma2 = sigma1; sigma3 = sigma1;

means    = zeros(1, 9);
std_devs = [sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3]; 
num_realizations = length(t); % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

% Integrate trajectory
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
STM0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;STM0], options);
% % rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;                                % constant position error
rn = state_t(:, 1:3)' + [sin(1E-4.*t);sin(1E-3.*t);sin(5E-4.*t)].*Ar;       % sinusoidal position error
vn = state_t(:, 4:6)';

% generate measurements
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm);

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma_n = [1E-2;1E-2;1E-2;1E-2;1E-2];
% % sigma_n = ones(10, 1).*1E-2;
% % [S] = mat2list(sigma_Cnm, sigma_Snm, Nc, Ns);
% % sigma_RMS     = computeRMS_coeffErr(n_max, Nc, Ns, S, Cnm.*0, Snm.*0);
% % sigma_n = sigma_RMS(2:end);

[Xp, Pp] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);
P0 = Pp(2:end, 2:end); 

% Gravity estimation
R_R = diag([sqrt(2)*sigma1, sigma2].^2);
R_N = diag([sigma1, sigma2, sigma3, sigma1, sigma2].^2);
n = ones(2, Nt)*NaN;

% loop
iterMax = 7;
count   = 0;
xnot_R = zeros(Ncs-1, 1); xnot_N = xnot_R; xnot_LS = xnot_R; xnot_N_AP = xnot_R;
Cp_R = Cp; Cp_N = Cp; Cp_LS = Cp; Cp_N_AP = Cp;
Sp_R = Sp; Sp_N = Sp; Sp_LS = Sp; Sp_N_AP = Sp;
Xp_R = Xp; Xp_N = Xp; Xp_LS = Xp; Xp_N_AP = Xp;

Pc = zeros(3, 3); Pxc = zeros(Ncs - 1, 3); Pxc_N_AP = zeros(Ncs - 1, 6);
Pc_N_AP = zeros(6, 6);
for j = 1:length(Pc)
    Pc(j, j) = Ar(j)^2;
end
for j = 1:length(Pc_N_AP)
    Pc_N_AP(j, j) = Ar(1)^4;
end
c = Ar.*1; % apriori values for the Consider Parameters; 
c_N_AP = ones(6, 1).*Ar(1)^2.*1;

while count < iterMax
    Ax_R = inv(P0); Ax_N = inv(P0); Ax = inv(P0); Ax_N_AP = inv(P0);
    Nx_R = -inv(P0) * xnot_R; Nx_N = -inv(P0) * xnot_N; Nx = -inv(P0) * xnot_LS;
    Nx_N_AP = -inv(P0) * xnot_N_AP;

    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    [~, Mxc_N_AP, Mcc_N_AP] = get_considerCov_apriori(P0, Pc_N_AP, Pxc_N_AP);
    
    for j = 1:Nt
        % RTN rotation matrix
        ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
        rn_RTN = ACI_RTN' * rn(:, j);
        rn_ACI = rn(:, j);
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % Rummel's method
        [Yc, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cp_R, Sp_R);

        [ax, nx] = rummels_method(Y(:, j)-Yc, Hc_RTN, R_R, ACI_RTN, noise(:, j));
        Ax_R  = Ax_R + ax;
        Nx_R  = Nx_R + nx;
    
        % Null space method
        [Yc, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cp_N, Sp_N);

        [Hpos] = compute_posPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);

        [ax, nx] = nullSpace_method(Y(:, j)-Yc, Hc_RTN, R_N, Hpos, ACI_RTN, noise(:, j));
        Ax_N  = Ax_N + ax;
        Nx_N  = Nx_N + nx;

        % Null space merthod + Apriori (AP)
        [Yc, Hc_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cp_N_AP, Sp_N_AP);

        [Hpos] = compute_posPartials(n_max, normalized, Cp_N_AP, Sp_N_AP, Re, GM, rn_ACI, ACAF_ACI);

        Hap = compute_posPartials_2ndOrder(GM, rn_ACI(1), rn_ACI(2), rn_ACI(3));
        [ax, nx, mxc, mcc] = nullSpace_method_AP(Y(:, j)-Yc, Hc_ACI, R_N, Hpos, Hap, eye(3,3), noise(:, j));
        Ax_N_AP  = Ax_N_AP + ax;
        Nx_N_AP  = Nx_N_AP + nx;
        Mxc_N_AP = Mxc_N_AP + mxc;
        Mcc_N_AP = Mcc_N_AP + mcc;
        
        % LS method
        [Yc, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cp_LS, Sp_LS);

        [Hpos] = compute_posPartials(n_max, normalized, Cp_LS, Sp_LS, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);

        [ax, nx, mxc, mcc] = LS_method(Y(:, j)-Yc, Hc_RTN, R_N, Hpos, ACI_RTN, noise(:, j));
        Ax = Ax + ax;
        Nx = Nx + nx;
        Mxc = Mxc + mxc;
        Mcc = Mcc + mcc;
    end

    % solve LS
    XNOT_R = Ax_R\Nx_R;
    XNOT_N = Ax_N\Nx_N;
    XNOT_LS = Ax\Nx - Ax\(Mxc * c);
    XNOT_N_AP = Ax_N_AP\Nx_N_AP - Ax_N_AP\(Mxc_N_AP * c_N_AP);

    Xp_R(2:end) = Xp_R(2:end) + XNOT_R;
    Xp_N(2:end) = Xp_N(2:end) + XNOT_N;
    Xp_LS(2:end) = Xp_LS(2:end) + XNOT_LS;
    Xp_N_AP(2:end) = Xp_N_AP(2:end) + XNOT_N_AP;

    [Cp_R, Sp_R] = list2mat(n_max, Nc, Ns, Xp_R);
    [Cp_N, Sp_N] = list2mat(n_max, Nc, Ns, Xp_N);
    [Cp_LS, Sp_LS] = list2mat(n_max, Nc, Ns, Xp_LS);
    [Cp_N_AP, Sp_N_AP] = list2mat(n_max, Nc, Ns, Xp_N_AP);

    % update corrections
    xnot_R = xnot_R + XNOT_R;
    xnot_N = xnot_N + XNOT_N;
    xnot_LS = xnot_LS + XNOT_LS;
    xnot_N_AP = xnot_N_AP + XNOT_N_AP;

    % show error
    disp('Rummels update = '    + string(vecnorm(XNOT_R)));
    disp('Null space update = ' + string(vecnorm(XNOT_N)));
    disp('LS update = ' + string(vecnorm(XNOT_LS)));
    disp('Null space + AP update = ' + string(vecnorm(XNOT_N_AP)));

    % update counter
    count = count + 1;
end
Px = inv(Ax);
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';
Pxc = Sxc * Pc;

Px_N_AP = inv(Ax_N_AP);
Sxc_N_AP = -Px_N_AP * Mxc_N_AP;
Pxx_N_AP = Px_N_AP + Sxc_N_AP*Pc_N_AP*Sxc_N_AP';
Pxc_N_AP = Sxc_N_AP * Pc_N_AP;

P_N =  inv(Ax_N);
P_R =  inv(Ax_R);
P_LS =  [Pxx, Pxc;Pxc', Pc];
P_N_AP =  [Pxx_N_AP, Pxc_N_AP;Pxc_N_AP', Pc_N_AP];
sigma_N = sqrt(diag(P_N));
sigma_R = sqrt(diag(P_R));
sigma_LS = sqrt(diag(P_LS));
sigma_N_AP = sqrt(diag(P_N_AP));

[Xp_R] = mat2list(Cp_R, Sp_R, Nc, Ns);
[Xp_N] = mat2list(Cp_N, Sp_N, Nc, Ns);
[Xp_LS] = mat2list(Cp_LS, Sp_LS, Nc, Ns);
[Xp_N_AP] = mat2list(Cp_N_AP, Sp_N_AP, Nc, Ns);

SH_R = Xp_R(2:end);
SH_N = Xp_N(2:end);
SH_LS = Xp_LS(2:end);
SH_N_AP = Xp_N_AP(2:end);

% plot trajectory
tt = 'Orbit radius along trajectory. T = ' + string(T./3600) + ' h';
plot_orbit(state_t, name, t./T ,Re, tt)

% % % plot SH estimation
% % tt1 = 'Estimation value. Cnm coefficients';
% % tt2 = 'Estimation value. Snm coefficients';
% % lgn =  {'truth','PMTM', 'NSM'};
% % ls  = '-'; mk = 'square'; 
% % plot_gravField(X, SH_R, SH_N, n_max, tt1, tt2, ls, mk, lgn);
% % 
% % % plot SH estimation
% % tt1 = 'Estimation value. Cnm coefficients';
% % tt2 = 'Estimation value. Snm coefficients';
% % ls  = '-'; mk = 'square'; 
% % lgn =  {'truth','LS', 'NSM'};
% % plot_gravField(X, SH_LS, SH_N, n_max, tt1, tt2, ls, mk, lgn);
% % 
% % % plot uncertainty
% % tt1 = 'Uncertainty SH. Cnm coefficients';
% % tt2 = 'Uncertainty SH. Snm coefficients';
% % lgn =  {'truth','PMTM', 'NSM'};
% % ls  = '-'; mk = 'square'; 
% % plot_gravField(X, 3.*sigma_R, 3.*sigma_N, n_max, tt1, tt2, ls, mk, lgn);
% % 
% % % plot uncertainty
% % tt1 = 'Uncertainty SH. Cnm coefficients';
% % tt2 = 'Uncertainty SH. Snm coefficients';
% % lgn =  {'truth','LS', 'NSM'};
% % ls  = '-'; mk = 'square'; 
% % plot_gravField(X, 3.*sigma_LS,3.*sigma_N, n_max, tt1, tt2, ls, mk, lgn);

% ploting difference
tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
lgn =  {'truth','PMTM', 'NSM', '3\sigma PMTM', '3\sigma NSM'};
mk = '*'; c1 = 'm'; c2 = 'b';
plot_error(X.*NaN, X(2:end) - SH_R, X(2:end) - SH_N, 3.*sigma_R, 3.*sigma_N, n_max, tt1, tt2, mk, lgn, c1, c2);

tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
lgn =  {'truth','LS', 'NSM', '3\sigma LS', '3\sigma NSM'};
c1 = 'g'; c2 = 'b';
plot_error(X, X(2:end) - SH_LS, X(2:end) - SH_N, 3.*sigma_LS, 3.*sigma_N, n_max, tt1, tt2, mk, lgn, c1, c2);

tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
lgn =  {'truth','NSM + AP', 'NSM', '3\sigma NSM + AP', '3\sigma NSM'};
c1 = "#57BDFF"; c2 = 'b';
plot_error(X, X(2:end) - SH_N_AP, X(2:end) - SH_N, 3.*sigma_N_AP, 3.*sigma_N, n_max, tt1, tt2, mk, lgn, c1, c2);

% % % plot correlation NSM
% % R_NSM = corrcov(P_N);
% % figure()
% % heatmap(R_NSM);


%% FUNCTIONS
function [ax, nx] = rummels_method(Y, Hc, R, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1), Y(2), Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;
    
    % select measurements
    dY = [dy(5) - dy(9); dy(8)];
    hc = [Hc(5, 2:end) - Hc(9, 2:end); Hc(8, 2:end)];
    
    % information and normal matrices
    ax = hc' * inv(R) * hc;
    nx = hc' * inv(R) * dY;
end

function [ax, nx] = nullSpace_method(Y, Hc, R, Hp, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    r  = C' * R * C;

% %     % de-correlate measurements
% %     [v, ~] = eig(r);
% %     r = v'*r*v;
% %     y = v'*y;
% %     hc = v'*hc;

    % information and normal matrices
    ax = hc' * inv(r) * hc;
    nx = hc' * inv(r) * y;
end

function [ax, nx, mxc, mcc] = nullSpace_method_AP(Y, Hc, R, Hp, Hap, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    hap = C' * Hap;
    r  = C' * R * C;

    % information and normal matrices
    ax = hc' * inv(r) * hc;
    nx = hc' * inv(r) * y;
    mxc = (hc' * inv(r) * hap);
    mcc = (hap' * inv(r) * hap);
end

function [ax, nx, mxc, mcc] = LS_method(Y, Hc, R, Hp, ACI_RTN, noise)
        % reshape meas [ACI]
        ddU = [Y(1), Y(2), Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
        % Rotate meas to RTN coordinates
        dy_RTN = ACI_RTN' * ddU * ACI_RTN;
        dy = reshape(dy_RTN, [9, 1]) + noise;

        hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
                Hc(8, 2:end)];
       
        hp = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)];
    
        % select measurements
        dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];
 
        ax  = hc' * inv(R) * hc;
        nx  = hc' * inv(R) * (dY);
        mxc = (hc' * inv(R) * hp);
        mcc = (hp' * inv(R) * hp); 
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


function [] = plot_orbit(state_t, name, time, Re, tt)
    plot_trajectory(state_t, name);
    Nt = length(time);
    figure()
    plot(time, vecnorm(state_t(:, 1:3)'), 'LineWidth', 2)
    hold on;
    plot(time, ones(1, Nt)*Re, 'LineWidth', 2, 'Color', 'r', 'LineStyle','--')
    xlabel('Orb. Period number, T')
    ylabel('[m]')
    title(tt)
    legend('orbit radius', 'brillouin sphere')
end

function [] = plot_gravField(X, SH_R, SH_N, n_max, tt1, tt2, ls, mk, lgn)
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    figure()
    subplot(1, 2, 1)
    semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Nc-1, abs(SH_R(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    title(tt1)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    
    subplot(1, 2, 2)
    semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold on;
    semilogy(1:Ns, abs(SH_R(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
    title(tt2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    legend(lgn);
end

function [] = plot_error(X, SH_R_err, SH_N_err, SH_R_sig, SH_N_sig, n_max, tt1, tt2, mk, lgn, c1, c2)
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    figure()
    subplot(1, 2, 1)
    semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Nc-1, abs(SH_R_err(1:Nc-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', c1, 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N_err(1:Nc-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', c2, 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_R_sig(1:Nc-1)), 'Marker','sq', 'LineStyle','-', 'LineWidth', 2, 'Color', c1, 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N_sig(1:Nc-1)), 'Marker','sq', 'LineStyle','-', 'LineWidth', 2, 'Color', c2, 'MarkerFaceColor', 'auto')
    title(tt1)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    
    subplot(1, 2, 2)
    semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Ns, abs(SH_R_err(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', c1, 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N_err(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle','--', 'LineWidth', 2, 'Color', c2, 'MarkerFaceColor','auto')
    semilogy(1:Ns, abs(SH_R_sig(Nc:Nc+Ns-1)), 'Marker','sq', 'LineStyle','-', 'LineWidth', 2, 'Color', c1, 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N_sig(Nc:Nc+Ns-1)), 'Marker','sq', 'LineStyle','-', 'LineWidth', 2, 'Color', c2, 'MarkerFaceColor','auto')
    title(tt2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    legend(lgn);
end