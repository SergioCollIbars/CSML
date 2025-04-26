clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Null space approach including attitude + position errors.
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

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.35E3;
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
f = 1/10;
t = linspace(0, rev*T, rev*T * f);
dt = t(2) - t(1);
Nt = length(t);

% position error
Ar = 1E-5*[1;1;1];                                 % [ACI]

% attitude error
Amp = 2E-9;
At      = Amp*t.*ones(3, Nt).*[1;0.7;0.5];
dA_dt   = Amp.*ones(3, Nt).*[1;0.7;0.5];        % [rad/s]
ddA_ddt = zeros(3, Nt);                         % [rad/s^2]

% attitude nominal value
Amp = 0;                                    % [rad]
attitude  = Amp.*t.*ones(3, Nt).*[0;1;0];   % nominal attitude [rad]
datt_dt   = Amp.*ones(3, Nt).*[0;1;0];
ddatt_ddt = zeros(3, Nt); 
[angVel_true, angAcc_true] = compute_angularVals(attitude + At, datt_dt + dA_dt, ddatt_ddt + ddA_ddt);
[angVel_nom, angAcc_nom]   = compute_angularVals(attitude, datt_dt, ddatt_ddt);

% plot S/C attitude
scale = 3600 * 180 / pi;
plot_Attitude(t./T, attitude, At * scale, angVel_nom, ...
    angVel_true - angVel_nom, angAcc_nom, angAcc_true - angAcc_nom);

% noise values from GOCE mission
noise0 = zeros(9, Nt);
% % sigma1  = 0.01 * 1E-9 * sqrt(f); % Vxx, Vyy
% % sigma2  = 0.6  * 1E-9 * sqrt(f); % Vyz, Vyx
% % sigma3  = 0.02 * 1E-9 * sqrt(f); % Vxz, Vzz
sigma1 = 1E-12;
sigma2 = sigma1; sigma3 = sigma1;

means    = zeros(1, 9);
std_devs = [sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3]; 
num_realizations = length(t); % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
rt = state_t(:, 1:3)';
vt = state_t(:, 4:6)';

% generate measurements & add rotations value
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm, eye(3,3));
[Y] = add_angularComponents(Y, attitude, At, flipud(angVel_true),...
    flipud(angAcc_true));  

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma_n = 1E-3*[1E-4;1E-4;1E-4;1E-2;1E-2];
% % sigma_n = ones(10, 1).*1E-4;
[Xp, ~] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_n] = ode113(@(t, x) EoM(t, x, Cp, Sp, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0 + Ar;v0;PHI0], options);
rn = state_n(:, 1:3)';
vn = state_n(:, 4:6)';

% plot position error
errorP = rt - rn;
figure()
plot(t./T, errorP, 'LineWidth', 2);
title('Position error over time')
ylabel('[m]')
xlabel('Orbital rev.')
legend('\Delta x', '\Delta y', '\Delta z')

% intial uncertainty
sigma_n = [1E0;1E0;1E0;1E0;1E0];
[~, Pp] = perturb_coeff(sigma_n, n_max, X);
P0 = Pp(2:end, 2:end); 

% Gravity estimation
R_NP = diag([sigma1, sigma2, sigma3, sigma1, sigma2].^2);
R_N = diag([sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3].^2);

Pc = zeros(3, 3); Pxc = zeros(Ncs - 1, 3);
for j = 1:3
    Pc(j, j) = 1.3*rms(errorP(j, :))^2;
end

c = sqrt(diag(Pc)).*1; % apriori values for the Consider Parameters;

% loop
iterMax = 10;
count   = 0;
iter    = 1;
xnot_NP = zeros(Ncs-1, 1); xnot_N = xnot_NP; xnot_LS = xnot_NP;
Cp_NP = Cp; Cp_N = Cp; Cp_LS = Cp;
Sp_NP = Sp; Sp_N = Sp; Sp_LS = Sp;
Xp_NP = Xp; Xp_N = Xp; Xp_LS = Xp;
rn_LS = rn; rn_N = rn; state_n_N = state_n;
while count < iterMax
    Ax_NP = inv(P0); Ax_N = inv(P0); Ax_LS = inv(P0);
    Nx_NP = -inv(P0) * xnot_NP; Nx_N = -inv(P0) * xnot_N; Nx_LS = -inv(P0) * xnot_LS;
    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    for j = 1:Nt
        % Position vector
        rn_ACI_N = rn_N(:, j);
        rn_ACI_LS = rn_LS(:, j);
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

        % from ACI to Nominal body frame
        B_ACI =rotationMatrix(attitude(1, j), attitude(2, j), attitude(3, j), ...
          [1, 2, 3]);   
        ACAF_B = ACAF_ACI * B_ACI';
    
         % compute attitude partials. Nominal body frame
        [Hrot_grad] = compute_rotPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_ACI_N, ACAF_ACI, ACAF_B);
        [Hrot_dA_ang, Hrot_dAdT_ang] = compute_angularDyadPartials(flipud(angVel_nom(:, j)), attitude(:, j), datt_dt(:, j), ddatt_ddt(:, j));
        Hrot_dA = Hrot_grad + Hrot_dA_ang;
        Hrot_dAdT = Hrot_dAdT_ang;
        [Hpos] = compute_posPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_ACI_N, ACAF_ACI, ACAF_B);
    
        % Null space method correcting for attitude
        [Yc, ~, Hc_BODY] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ACI_N', zeros(1, 3)], ...
                noise0, Cp_N, Sp_N, B_ACI');
        [Yc] = add_angularComponents(Yc, attitude(:, j), zeros(3, Nt), flipud(angVel_nom(:, j)),...
            flipud(angAcc_nom(:, j)));
      
        if(iter == 1)
            PHIi0 = reshape(state_n_N(j, 7:end), [6,6]);
            Hpr = [Hpos, Hrot_dA, Hrot_dAdT];
            dY = Y(:, j) - Yc + noise(:, j);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
                    Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];

            % update iter
            iter = 2;
        else
            PHIj0 = reshape(state_n_N(j, 7:end), [6,6]);
            PHIji = PHIj0/(PHIi0);
            PHI = PHIji(1:3, 1:3);
            Hrot = [Hpos, Hrot_dA, Hrot_dAdT]*[PHI, zeros(3, 6);...
                zeros(3, 3), eye(3,3),0.*eye(3,3);...
                zeros(3, 3), zeros(3,3), eye(3,3)];

            dy =  Y(:, j) - Yc + noise(:, j);
            dY = [dY; dy];

            hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(2, 2:end); Hc_BODY(5, 2:end);...
                    Hc_BODY(8, 2:end);  Hc_BODY(3, 2:end); Hc_BODY(6, 2:end); Hc_BODY(9, 2:end)];
            Hc = [Hc; hc];

            hpr = Hrot;
            Hpr = [Hpr; hpr];

            C   = null(Hpr');
            R = blkdiag(R_N, R_N);
            
            % project into null space
            r  = C' * R * C;
            y  = C' * dY;
            hc = C' * Hc;

            % information and normal matrices
            ax = hc' * inv(r) * hc;
            nx = hc' * inv(r) * y;

            Ax_N  = Ax_N + ax;
            Nx_N  = Nx_N + nx;

            % re-set iter
            iter = 1;
        end

        % classic LS
        [Yc, ~, Hc_BODY] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ACI_LS', zeros(1, 3)], ...
                noise0, Cp_LS, Sp_LS, B_ACI');
        [Yc] = add_angularComponents(Yc, attitude(:, j), zeros(3, 1), flipud(angVel_nom(:, j)),...
            flipud(angAcc_nom(:, j)));

        [Hpos] = compute_posPartials(n_max, normalized, Cp_LS, Sp_LS, Re, GM, rn_ACI_LS, ACAF_ACI, ACAF_B);

        [ax, nx, mxc, mcc] = LS_method(Y(:, j), Yc, Hc_BODY, Hpos, R_N, noise(:, j));
        Ax_LS  = Ax_LS + ax;
        Nx_LS  = Nx_LS + nx;
        Mxc = Mxc + mxc;
        Mcc = Mcc + mcc;

    end

    % solve LS
    XNOT_N = Ax_N\Nx_N;
    XNOT_LS = Ax_LS\Nx_LS;

    Xp_N(2:end) = Xp_N(2:end) + XNOT_N;
    Xp_LS(2:end) = Xp_LS(2:end) + XNOT_LS;

    [Cp_N, Sp_N] = list2mat(n_max, Nc, Ns, Xp_N);
    [Cp_LS, Sp_LS] = list2mat(n_max, Nc, Ns, Xp_LS);

    % update corrections
    xnot_N = xnot_N + XNOT_N;
    xnot_LS = xnot_LS + XNOT_LS;

    % show error
    disp('Null space update, attitude + position = ' + string(vecnorm(XNOT_N)));
    disp('Least Squares update = ' + string(vecnorm(XNOT_LS)));
    disp('--------------------------------------------------------'); 

    % recompute STM for the given orbit
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    PHI0 = reshape(eye(6,6), [36, 1]);
    [~, state_n_N] = ode113(@(t, x) EoM(t, x, Cp_N, Sp_N, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0 + Ar;v0;PHI0], options);
    rn_N = state_n_N(:, 1:3)';

% %     [STM] = integrateSTM(t, state_n, asterParams, poleParams, Cp_N, Sp_N, 5);
% %     state_n_N = [zeros(Nt, 6), STM];

    % update counter
    count = count + 1;
end
Px = inv(Ax_LS);
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';
Pxc = Sxc * Pc;

P_N =  inv(Ax_N);
P_LS =  [Pxx, Pxc;Pxc', Pc];
sigma_N = sqrt(diag(P_N));
sigma_LS = sqrt(diag(P_LS));

[Xp_N] = mat2list(Cp_N, Sp_N, Nc, Ns);
[Xp_LS] = mat2list(Cp_LS, Sp_LS, Nc, Ns);

SH_N = Xp_N(2:end);
SH_LS = Xp_LS(2:end);


% plot trajectory
tt = 'Orbit radius along trajectory. T = ' + string(T./3600) + ' h';
plot_orbit(state_t, name, t./T ,Re, tt)


% plot uncertainty
tt1 = 'Uncertainty SH. Cnm coefficients';
tt2 = 'Uncertainty SH. Snm coefficients';
ls  = '-'; mk = 'square'; 
llg = {'truth','LS', 'NSM'};
plot_gravField(X, 3.*sigma_LS(1:45), 3.*sigma_N, n_max, tt1, tt2, ls, mk, llg);


tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
ls  = '--'; mk = '*'; 
llg = {'truth','LS', 'NSM'};
plot_gravField(X.*NaN, X(2:end) - SH_LS, X(2:end) - SH_N, n_max, tt1, tt2, ls, mk, llg);

%% FUNCTIONS

function [ax, nx, mxc, mcc] = LS_method(Y, Yc, Hc, Hpos, R, noise)
    % add noise
    Y = Y + noise(1:9);

% %     % select measurements
% %     dY = [Y(1)-Yc(1);Y(2)-Yc(2);Y(3)-Yc(3);Y(5)-Yc(5);Y(6)-Yc(6)];
% % 
% % 
% %     % project measurements
% %     hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
% %         Hc(8, 2:end)];
    dY = Y - Yc;
    hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end); Hc(2, 2:end); Hc(5, 2:end);...
             Hc(8, 2:end);  Hc(3, 2:end); Hc(6, 2:end); Hc(9, 2:end)];
    hrp = Hpos;

    % information and normal matrices
    ax = hc' * inv(R) * hc;
    nx = hc' * inv(R) * dY;
    mxc = (hc' * inv(R) * hrp);
    mcc = (hrp' * inv(R) * hrp); 
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

function [] = plot_gravField(X, SH_R, SH_N, n_max, tt1, tt2, ls, mk, llg)
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
    legend(llg)
end

function [] = plot_Attitude(t, attitude_nom, attitude_true, angVel_nom, ...
    angVel_true, angAcc_nom, angAcc_true)
    figure()
    subplot(3, 2, 1)
    plot(t, attitude_true, 'LineWidth', 2)
    legend('\delta \Psi, yaw', '\delta \theta, pitch', '\delta \phi roll')
    ylabel('Arcsecond')
    title('Error attitude values')
    subplot(3, 2, 3)
    plot(t, rad2deg(angVel_true), 'LineWidth', 2)
    ylabel('[deg/s]')
    legend({'$\delta \omega_3$', '$\delta \omega_2$', '$\delta \omega_1$'}, 'Interpreter', 'latex')
    subplot(3, 2, 5)
    plot(t, rad2deg(angAcc_true), 'LineWidth', 2)
    ylabel('[deg/s^2]')
    legend({'$\delta \dot{\omega_3}$', '$\delta \dot{\omega_2}$', '$\delta \dot{\omega_1}$'}, 'Interpreter', 'latex')
    xlabel('Orbit revolution')
    subplot(3, 2, 2)
    f = wrapTo2Pi(attitude_nom);
    plot(t, rad2deg(f), 'LineWidth', 2)
    legend('\Psi, yaw', '\theta, pitch', '\phi roll')
    ylabel('[deg]')
    title('Nominal attitude values')
    subplot(3, 2, 4)
    plot(t, rad2deg(angVel_nom), 'LineWidth', 2)
    ylabel('[deg/s]')
    legend({'$\omega_3$', '$\omega_2$', '$\omega_1$'}, 'Interpreter', 'latex')
    subplot(3, 2, 6)
    plot(t, rad2deg(angAcc_nom), 'LineWidth', 2)
    ylabel('[deg/s^2]')
    legend({'$\dot{\omega_3}$', '$\dot{\omega_2}$', '$\dot{\omega_1}$'}, 'Interpreter', 'latex')
    xlabel('Orbit revolution')
end