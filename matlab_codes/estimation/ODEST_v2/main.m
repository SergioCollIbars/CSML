clear;
clc;
close all;
format long g;

set(0,'defaultAxesFontSize',16);

addpath('data/');
addpath('functions/integrator/');
addpath('functions/measurements/');
addpath('functions/plotter/');
addpath('functions/solvers/');
addpath('../../../QGG_gravEstim/data_files/') ;
addpath('../../../QGG_gravEstim/src/')

%%                         ODEST ALGORITHM V2
% Description: Improved version of the ODEST algorithm based on Rummel's
% method. Look for null space of sensitivity matrices to decouple position
% and gravity field errors in measurement space.


% Simulation parameters
n_max = 6;                 % max SH degree
normalized = 1;            % Normalized coefficients
sigma = 6.3E-12;           % Noise level
f = 1/10;                  % Meas frequency [Hz]

% Harmonics values. Truth (matrix form)
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm_t, Snm_t, Re] = readCoeff(path);

% Harmonic values. A priori (vector form)
[Nc, Ns, Ncs] = count_num_coeff(n_max);
[Cnm_ap, Snm_ap] = getkaula(n_max, Nc, Ns, normalized);
[Xt] = mat2list(Cnm_t, Snm_t, Nc, Ns);

% Asteroid parameters. Bennu
GM =  5.2;                % Point mass value    [m^3/s^2]
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % % Asteroid parameters. Eros
% % GM =  459604.431484721;          % Point mass value    [m^3/s^2]
% % W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
% % W0 = 0;                          % Initial asteroid longitude
% % RA = deg2rad(11.363);            % Right Ascension     [rad]
% % DEC = deg2rad(17.232);           % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% Initial trajectory & time points
r = 0.3E3;                                        % orbit radius    [m]
phi = pi/2;                                     % orbit latitude  [rad]
lambda = 0;                                  % orbit longitude [rad]
theta = pi/2 - phi;                             % orbit colatitude [rad]

R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];
v0 = R * [0;0;sqrt(GM/r)];

% % r0 = [3.6998; -105.0537; 994.459667937083];     % Position            [m]
% % v0 = [-0.05286612;-0.04879135;-0.004957585];    % Velocity            [m]

% % r0 = [ 7.39956897047445;-210.107405307177;1988.91933587417];
% % v0 = [-0.0373819912839047;-0.034500694215296; -0.00350554222806066];

n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = 1 * (2 * pi / n);
t = linspace(0, T, T * f);
Nt = length(t);

% Real trajectory & plot
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm_t, Snm_t, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;reshape(eye(6,6), [36, 1])], options);

plot_trajectory(state_t, "BENNU");

% generate noise
% % noise = normrnd(0,sigma, [9, length(t)]);
noise = load('noise_6E12_r_1km.mat').noise;

% generate measurements
[ddU_ACI, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise, Cnm_t, Snm_t);

% nominal position
Ar = 0.5*[1;1;1];             % position error ACI [m]
Av = [5;5;9]*0;               % velocity error ACI [m/s]

rn0 = state_t(1, 1:3)' + Ar;
vn0 = state_t(1, 4:6)' + Av;

% ODEST V2 algorithm
iterMax = 10;
count = 0;

eps = 1E-4;
err = eps + 1;

% transform to matrix 
[C, S] = list2mat(n_max, Nc, Ns, [Cnm_ap;Snm_ap]);

% apriori uncertanties
sigmaC = 1;
sigmaS = 1;
sigmaP = 10;
sigmaV = 1E-2;
P0_c = diag([repelem(sigmaC, Nc-1), repelem(sigmaS, Ns)].^2);
P0_p = diag([sigmaP, sigmaP, sigmaP, sigmaV, sigmaV, sigmaV].^2);
xnot_c = zeros(Nc+Ns-1, 1);
xnot_p = zeros(6, 1);

% integrate nominal trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6,6), [36, 1]);
[~, state] = ode113(@(t, x) EoM(t, x, C, S, n_max, GM, Re, ...
    normalized, W0, W, RA, DEC), t, [rn0;vn0;STM0], options);

% nominal states
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;
vn = state_t(:, 4:6)';

PHI_t = state(:, 7:end);

% output values
errState = ones(6, Nt, iterMax) * NaN;
sigmaState = ones(6, Nt, iterMax) * NaN;
RMSerr   = ones(iterMax, n_max) * NaN;
RMSsigma = ones(iterMax, n_max) * NaN;
coeffErr = ones(iterMax, Nc+Ns) * NaN;

% init values
errState(:, :, 1) = [rn;vn] - state_t(:, 1:6)';
[X] = mat2list(C, S, Nc, Ns);
RMSerr(1, :) = computeRMS_coeffErr(n_max, Nc, Ns, ...
    X, Cnm_t, Snm_t);
coeffErr(1, :) = abs((X-Xt)./Xt);
st = zeros(n_max+1, n_max+1);
RMSsigma(1, :) = computeRMS_coeffErr(n_max, Nc, Ns, ...
    [0; sqrt(diag(P0_c))], st, st);
for j = 1:Nt
    PHI_i0 = reshape(PHI_t(j, :), [6,6]);
    s = PHI_i0 * P0_p * PHI_i0';
    sigmaState(:, j, 1) = sqrt(diag(s));
end

% loop
while (err> eps) && (count < iterMax)
    % run ODEST v2
    [Cnm_new, Snm_new, state_new, XNOT_c, XNOT_p, Pc_new, Pp_new, PHI] =...
        ODEST_v2(t ,ddU_ACI, rn, vn, C, S, poleParams, asterParams, ...
        sigma, xnot_c, xnot_p, P0_c, P0_p);

% %     % run ODEST v4
% %     [Cnm_new, Snm_new, state_new, XNOT_c, XNOT_p, Pc_new, Pp_new, PHI] =...
% %     ODEST_v4(t ,ddU_ACI, rn, vn, C, S, poleParams, asterParams, ...
% %     sigma, xnot_c, xnot_p, P0_c, P0_p);
% % 
% %     % run ODEST v0
% %     [Cnm_new, Snm_new, state_new] =ODEST_v0(t ,ddU_ACI, rn, vn, C, S, poleParams, asterParams, ...
% %     sigma, P0_c);

    % update deviations
    xnot_c = xnot_c + XNOT_c;
    xnot_p = xnot_p + XNOT_p;
    
    % update error
    err = vecnorm(XNOT_p(1:3));

    % update inital conditions
    rn = state_new(1:3, :);
    vn = state_new(4:6, :);
    C = Cnm_new;
    S = Snm_new;

    % update counter
    count = count + 1;

    % save state errors
    errState(:, :, count+1) = [rn;vn] - state_t(:, 1:6)';
    for j = 1:Nt
        PHI_i0 = reshape(PHI(j, :), [6,6]);
        s = PHI_i0 * Pp_new * PHI_i0';
        sigmaState(:, j, count+1) = sqrt(diag(s));
    end

    % compute RMS coefficient error
    [X] = mat2list(Cnm_new, Snm_new, Nc, Ns);
    RMSerr(count+1, :) = computeRMS_coeffErr(n_max, Nc, Ns, ...
        X, Cnm_t, Snm_t);
    RMSsigma(count+1, :) = computeRMS_coeffErr(n_max, Nc, Ns, ...
        [0; sqrt(diag(Pc_new))], st, st);
    coeffErr(count+1, :) = abs((X-Xt)./Xt);
    disp(count)
end

% plot SH errors
plot_coeffErr(n_max, RMSerr, RMSsigma, coeffErr,count);

% plot state errors
plot_stateErr(errState, sigmaState, t, count)


[X] = mat2list(Cnm_t, Snm_t, Nc, Ns);
SH_R = mat2list(Cnm_new, Snm_new, Nc, Ns); 
[num_C, num_S, str_C, str_S] = SH_xlabel(n_max);

figure()
subplot(1, 2, 1)
semilogy(1:Nc-1, abs(X(1:Nc-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold all;
semilogy(1:Nc-1, abs(SH_R(1:Nc-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
title('Estimation SH Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;

subplot(1, 2, 2)
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold on;
semilogy(1:Ns, abs(SH_R(Nc:Nc+Ns-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
title('Estimation SH. Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;
legend('truth','n=0 truncation', 'null-space')


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