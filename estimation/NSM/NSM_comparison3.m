clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Null space approach including attitude errors.
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

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.3E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 0.5*[1;1;1];            % [ACI]

% attitude error
At = 5E-5.*[1;1;1];         % [yaw, pitch, roll]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% noise values from GOCE mission
noise0 = zeros(9, Nt);
sigma1  = 0.01 * 1E-9 * sqrt(f); % Vxx, Vyy
sigma2  = 0.6  * 1E-9 * sqrt(f); % Vyz, Vyx
sigma3  = 0.02 * 1E-9 * sqrt(f); % Vxz, Vzz
sigmaq  = 5E-6;                  % rad 
% % sigma2 = sigma1; sigma3 = sigma1;

means    = zeros(1, 12);
std_devs = [sigma1, sigma2, sigma3, sigma2, sigma1, sigma2, sigma3, ...
    sigma2, sigma3, sigmaq, sigmaq, sigmaq]; 
num_realizations = length(t); % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;PHI0], options);
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;
vn = state_t(:, 4:6)';

% generate measurements
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm);

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma_n = [1E-2;1E-2;1E-2;1E-2;1E-2];
% % sigma_n = ones(10, 1).*1E-4;
[Xp, Pp] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);
P0 = Pp(2:end, 2:end); 

% Gravity estimation
R_N = diag([sigma1, sigma2, sigma3, sigma1, sigma2, sigmaq, sigmaq, sigmaq].^2);
R_NP = diag([sigma1, sigma2, sigma3, sigma1, sigma2].^2);
n = ones(2, Nt)*NaN;

% loop
iterMax = 10;
count   = 0;
xnot_NP = zeros(Ncs-1, 1); xnot_N = xnot_NP;
Cp_NP = Cp; Cp_N = Cp;
Sp_NP = Sp; Sp_N = Sp;
Xp_NP = Xp; Xp_N = Xp;
while count < iterMax
    Ax_NP = inv(P0); Ax_N = inv(P0);
    Nx_NP = -inv(P0) * xnot_NP; Nx_N = -inv(P0) * xnot_N;
    for j = 1:Nt
        % RTN rotation matrix
        ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
        rn_RTN = ACI_RTN' * rn(:, j);
        rn_ACI = rn(:, j);
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
         % compute position and attitude partials
        [Hpos] = compute_posPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);

        [Hrot] = compute_rotPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);
    
        % Null space method correcting for attitude
        [Yc, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
                noise0, Cp_N, Sp_N);

        [ax, nx] = nullSpace_method2(Y(:, j), Yc, [Hc_RTN;zeros(3, Ncs)], R_N, [Hpos,Hrot], ACI_RTN, noise(:, j), At);
        Ax_N  = Ax_N + ax;
        Nx_N  = Nx_N + nx;

        % Null space method only accounting for position
        [Hpos] = compute_posPartials(n_max, normalized, Cp_NP, Sp_NP, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);
        [ax, nx] = nullSpace_method(Y(:, j), Yc, Hc_RTN, R_NP, Hpos, ACI_RTN, noise(:, j), At);
        Ax_NP  = Ax_NP + ax;
        Nx_NP  = Nx_NP + nx;
    end

    % solve LS
    XNOT_NP = Ax_NP\Nx_NP;
    XNOT_N = Ax_N\Nx_N;

    Xp_NP(2:end) = Xp_NP(2:end) + XNOT_NP;
    Xp_N(2:end) = Xp_N(2:end) + XNOT_N;

    [Cp_NP, Sp_NP] = list2mat(n_max, Nc, Ns, Xp_NP);
    [Cp_N, Sp_N] = list2mat(n_max, Nc, Ns, Xp_N);

    % update corrections
    xnot_NP = xnot_NP + XNOT_NP;
    xnot_N = xnot_N + XNOT_N;

    % show error
    disp('Null space for position update = '  + string(vecnorm(XNOT_NP)));
    disp('Null space update = ' + string(vecnorm(XNOT_N)));
    
    % update counter
    count = count + 1;
end

P_N =  inv(Ax_N);
P_NP =  inv(Ax_NP);
sigma_N = sqrt(diag(P_N));
sigma_NP = sqrt(diag(P_NP));

[Xp_NP] = mat2list(Cp_NP, Sp_NP, Nc, Ns);
[Xp_N] = mat2list(Cp_N, Sp_N, Nc, Ns);

SH_NP = Xp_NP(2:end);
SH_N = Xp_N(2:end);


% plot trajectory
tt = 'Orbit radius along trajectory. T = ' + string(T./3600) + ' h';
plot_orbit(state_t, name, t./T ,Re, tt)

% plot SH estimation
tt1 = 'Estimation value. Cnm coefficients';
tt2 = 'Estimation value. Snm coefficients';
ls  = '-'; mk = 'square'; 
plot_gravField(X, SH_NP, SH_N, n_max, tt1, tt2, ls, mk);

% plot uncertainty
tt1 = 'Uncertainty SH. Cnm coefficients';
tt2 = 'Uncertainty SH. Snm coefficients';
ls  = '-'; mk = 'square'; 
plot_gravField(X, sigma_NP, sigma_N, n_max, tt1, tt2, ls, mk);

% ploting difference
tt1 = 'Estimation error. Cnm coefficients';
tt2 = 'Estimation error. Snm coefficients';
ls  = '--'; mk = '*'; 
plot_gravField(X.*NaN, X(2:end) - SH_NP, X(2:end) - SH_N, n_max, tt1, tt2, ls, mk);

%% FUNCTIONS
function [ax, nx] = nullSpace_method(Y, Yc, Hc, R, Hp, ACI_RTN, noise, At)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
    ddUc = [Yc(1),Yc(2),Yc(3);Yc(4),Yc(5),Yc(6);Yc(7),Yc(8),Yc(9)];
    [Rot] = rotationMatrix(At(1), At(2), At(3), [3, 2, 1]);

    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy_RTN = Rot' * dy_RTN * Rot;

    dy = reshape(dy_RTN, [9, 1]) + noise(1:9);
    dyc = reshape(ACI_RTN' * ddUc * ACI_RTN, [9, 1]);

    % select measurements
    dY = [dy(1)-dyc(1);dy(4)-dyc(4);dy(7)-dyc(7);dy(5)-dyc(5);dy(8)-dyc(8)];

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    r  = C' * R * C;

    % de-correlate measurements
    [v, ~] = eig(r);
    r = v'*r*v;
    y = v'*y;
    hc = v'*hc;

    % information and normal matrices
    ax = hc' * inv(r) * hc;
    nx = hc' * inv(r) * y;
end

% NSM accounting for position and attitude errors
function [ax, nx] = nullSpace_method2(Y, Yc, Hc, R, Hpr, ACI_RTN, noise, At)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
    ddUc = [Yc(1),Yc(2),Yc(3);Yc(4),Yc(5),Yc(6);Yc(7),Yc(8),Yc(9)];
    [Rot] = rotationMatrix(At(1), At(2), At(3), [3, 2, 1]);

    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy_RTN = Rot' * dy_RTN * Rot;

    dy = reshape(dy_RTN, [9, 1]) + noise(1:9);
    dyc = reshape(ACI_RTN' * ddUc * ACI_RTN, [9, 1]);

    % select measurements
    dY = [dy(1)-dyc(1);dy(4)-dyc(4);dy(7)-dyc(7);dy(5)-dyc(5);dy(8)-dyc(8);At+noise(10:12)];

    % look for null space
    H = [Hpr;zeros(3,3), eye(3,3)];
    C = null([H(1, :);H(2,:);H(3,:);H(5, :);H(6, :);H(10:12, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end); Hc(10:12, 2:end)];
    r  = C' * R * C;

    % de-correlate measurements
    [v, ~] = eig(r);
    r = v'*r*v;
    y = v'*y;
    hc = v'*hc;

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

function [] = plot_gravField(X, SH_R, SH_N, n_max, tt1, tt2, ls, mk)
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
    legend('truth','NSM', 'ENSM')
end