clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM VS THE LEAST SQUARES ALGORITHM VS EKF
% Description: Compare the NSM vs the LS and EKF estimating gravity field and
% position errors.
% Author: Sergio Coll
% Date: 10/22/2024

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

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 1E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 2E-3.*[1;1;1];            % [ACI]
Av = 0.*[1;1;1];           % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 5;
f = 1/30;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% noise values from GOCE mission
noise0 = zeros(9, Nt);
sigma1  = 0.01 * 1E-9 * sqrt(f); % Vxx, Vyy
sigma2  = 0.6  * 1E-9 * sqrt(f); % Vyz, Vyx
sigma3  = 0.02 * 1E-9 * sqrt(f); % Vxz, Vzz

sigma1 = 1E-13;
sigma2 = sigma1; sigma3 = sigma1;

means    = zeros(1, 9);
std_devs = [sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3]; 
num_realizations = length(t); % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma_n = [1E-6;1E-6;1E-6;1E-6;1E-6];
[Xp, ~] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;reshape(eye(6,6), [36, 1])], options);
[~, state_n] = ode113(@(t, x) EoM(t, x, Cp, Sp, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0+Ar;v0+Av;reshape(eye(6,6), [36, 1])], options);
rn = state_n(:, 1:3);
vn = state_n(:, 4:6);

figure()
plot(t, (state_t(:, 1:3)-state_n(:, 1:3))')
title('Position error in time');

% generate measurements
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm);

% Gravity estimation weight meas. Initial uncertainty
R_N = diag([sigma1, sigma2, sigma3, sigma1, sigma2].^2);

sigma_n0 = [1E-2;1E-2;1E-2;1E-2;1E-2];      % apriori uncertanty gra. field
sigmaPos = [1;1;1;1E-2;1E-2;1E-2];          % apriori uncertainty s/c state
[~, Pp] = perturb_coeff(sigma_n0, n_max, X);
P0 = Pp(2:end, 2:end); 

% define output values & estimation parameters
iterMax = 7;
count   = 0;
xnot_L = zeros(Ncs-1 + 6, 1); xnot_N = zeros(Ncs-1, 1);
Cp_L = Cp; Cp_N = Cp; Cp_E = Cp;
Sp_L = Sp; Sp_N = Sp; Sp_E = Sp;
Xp_N = Xp;
% solve for NSM
disp('Solve NSM')
while count < iterMax
    Ax_N = inv(P0);
    Nx_N = -inv(P0) * xnot_N;
    for j = 1:Nt
        % Position and velocity vector (used in NSM)
        rn_ACI = rn(j, :)';
        vn_ACI = vn(j, :)';

        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % Null space method
        [Hpos] = compute_posPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_ACI, ACAF_ACI);
        [Yc, HC_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ACI', vn_ACI'], ...
                noise0, Cp_N, Sp_N);

        [ax, nx] = nullSpace_method(Y(:, j)-Yc, HC_ACI, R_N, Hpos, eye(3,3), noise(:, j));
        Ax_N  = Ax_N + ax;
        Nx_N  = Nx_N + nx;
    end

    % solve LS with NSM
    XNOT_N = Ax_N\Nx_N;

    Xp_N(2:end) = Xp_N(2:end) + XNOT_N;

    [Cp_N, Sp_N] = list2mat(n_max, Nc, Ns, Xp_N);

    % update corrections
    xnot_N = xnot_N + XNOT_N;

    % show error
    disp('Null space update = ' + string(vecnorm(XNOT_N)));

    % update counter
    count = count + 1;
end


% solve for LS
disp('Solve LS')
[~, state_n] = ode113(@(t, x) EoM(t, x, Cp_L, Sp_L, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0+Ar;v0+Av;reshape(eye(6,6), [36, 1])], options);
STM = state_n(:, 7:end);
Xp_L = [Xp;state_n(1, 1:6)']; Xp_E = Xp_L;
count = 0;
iterMax = 20;
err = 0;
obs1 = ones(1, Nt) * NaN;
while (count < iterMax) && (err < 5)
    Ax_L = inv(blkdiag(P0, diag(sigmaPos.^2)));
    Nx_L = -inv(blkdiag(P0, diag(sigmaPos.^2))) * xnot_L;
    for j = 1:Nt        
        % STM at current time
        PHIt = reshape(STM(j, :), [6,6]);

        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % LS method
        [Hpos] = compute_posPartials(n_max, normalized, Cp_L, Sp_L, Re, GM, state_n(j, 1:3)', ACAF_ACI);
        [Yc, HC_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, state_n(j, 1:6), ...
                noise0, Cp_L, Sp_L);
        [ax, nx] = LS_method(Y(:, j)-Yc, HC_ACI, Hpos, R_N, PHIt, noise(:, j));

        Ax_L  = Ax_L + ax;
        Nx_L  = Nx_L + nx;
    end

    % solve LS
    XNOT_L = Ax_L\Nx_L;

    Xp_L(2:end) = Xp_L(2:end) + XNOT_L;

    [Cp_L, Sp_L] = list2mat(n_max, Nc, Ns, Xp_L(1:46));

    % update corrections
    xnot_L = xnot_L + XNOT_L;

    % compute new state
    r0 = Xp_L(47:49);
    v0 = Xp_L(50:52); 
    [~, state_n] = ode113(@(t, x) EoM(t, x, Cp_L, Sp_L, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;reshape(eye(6,6), [36, 1])], options);
    STM = state_n(:, 7:end);

    % show error
    err = vecnorm(XNOT_L);
    disp('LS update = '    + string(err));

    % update counter
    count = count + 1;
end

% solve for EKF
disp('Solve EKF')
P = blkdiag(P0, diag(sigmaPos.^2));
s = Xp_E(47:end);
for j = 2:Nt
    % time span
    t_span = [t(j - 1), t(j)];

    % integrate trajectory
     [~, state_n] = ode113(@(t, x) EoM(t, x, Cp_E, Sp_E, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t_span, [s;reshape(eye(6,6), [36, 1])], options);
    STM = state_n(:, 7:end);

    % current states 
    rn = state_n(end, 1:3);
    vn = state_n(end, 4:6);

    PHIt = reshape(STM(end, :), [6,6]);

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % LS method
    [Hpos] = compute_posPartials(n_max, normalized, Cp_E, Sp_E, Re, GM, rn', ACAF_ACI);
    [Yc, HC_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn, vn], ...
            noise0, Cp_E, Sp_E);
    [Xhat, P] = EKF_method(Y(:, j)-Yc, HC_ACI, Hpos, R_N, PHIt, noise(:, j), P);
    s = [rn,vn]' + Xhat(46:end);
    Xp_E(2:end) = Xp_E(2:end) + Xhat;
    [Cp_E, Sp_E] = list2mat(n_max, Nc, Ns, Xp_E(1:46));
end

P_N =  inv(Ax_N);
P_L =  inv(Ax_L);
P_E = P;
sigma_N = sqrt(diag(P_N));
sigma_L = sqrt(diag(P_L));
sigma_E = sqrt(diag(P_E));

[xp_L] = mat2list(Cp_L, Sp_L, Nc, Ns);
[xp_N] = mat2list(Cp_N, Sp_N, Nc, Ns);
[xp_E] = mat2list(Cp_E, Sp_E, Nc, Ns);

SH_L = Xp_L(2:46);
SH_N = Xp_N(2:end);
SH_E = Xp_E(2:46);

% plot trajectory
tt = 'Orbit radius along trajectory. T = ' + string(T./3600) + ' h';
plot_orbit(state_t, name, t./T ,Re, tt)

% plot SH estimation
tt1 = 'Estimation value. Cnm coefficients';
tt2 = 'Estimation value. Snm coefficients';
ls  = '-'; mk = 'square'; 
plot_gravField(X, SH_L, SH_N, SH_E, n_max, tt1, tt2, ls, mk);

% plot uncertainty
tt1 = 'Uncertainty SH. Cnm coefficients';
tt2 = 'Uncertainty SH. Snm coefficients';
ls  = '-'; mk = 'square'; 
plot_gravField(X, sigma_L, sigma_N, sigma_E, n_max, tt1, tt2, ls, mk);

% ploting difference
tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
ls  = '--'; mk = '*'; 
plot_gravField(X.*NaN, X(2:end) - SH_L, X(2:end) - SH_N,  X(2:end) - SH_E, n_max, tt1, tt2, ls, mk);


%% FUNCTIONS
function [X_hat, P] = EKF_method(Y, Hc, Hp, R0, PHI, noise, P)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy = reshape(ddU, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];
    
    % STM accounting for gravity field
    Nc  = length(Hc(1, 2:end));
    PHIt = [eye(Nc, Nc), zeros(Nc, 6); ...
        zeros(6, Nc), PHI];
    
    % compute measurement partials
    hp  = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)];
    hc  = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    ht  = [hc, hp, zeros(5, 3)] * PHIt;
    
    % information and normal matrices
    P_bar = PHIt * P * PHIt';
    K = P_bar * ht'/(ht * P_bar * ht' + R0);
    Nx = length(ht(1, :));

    % update meas
    X_hat = K * dY;
    P = (eye(Nx, Nx) - K * ht) * P_bar * (eye(Nx, Nx) - K * ht)' + ...
        K * R0 * K';
end

function [ax, nx] = LS_method(Y, Hc, Hp, R, PHI, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy = reshape(ddU, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];
    
    % compute measurement partials
    Nc  = length(Hc(1, 2:end));
    hp  = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)];
    hc  = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    ht  = [hc, hp, zeros(5, 3)] * [eye(Nc, Nc), zeros(Nc, 6); ...
        zeros(6, Nc), PHI];
    
    % information and normal matrices
    ax = ht' * inv(R) * ht;
    nx = ht' * inv(R) * dY;
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

function [] = plot_gravField(X, SH_R, SH_N, SH_E, n_max, tt1, tt2, ls, mk)
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
    figure()
    subplot(1, 2, 1)
    semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold all;
    semilogy(1:Nc-1, abs(SH_R(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_N(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
    semilogy(1:Nc-1, abs(SH_E(1:Nc-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor', 'auto')
    title(tt1)
    xticks(num_C);
    xticklabels(str_C);
    grid on;
    
    subplot(1, 2, 2)
    semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
    hold on;
    semilogy(1:Ns, abs(SH_R(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
    semilogy(1:Ns, abs(SH_N(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
    semilogy(1:Ns, abs(SH_E(Nc:Nc+Ns-1)), 'Marker',mk, 'LineStyle',ls, 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor','auto')
    title(tt2)
    xticks(num_S);
    xticklabels(str_S);
    grid on;
    legend('truth','LS', 'NSM', 'EKF')
end