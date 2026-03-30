clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Compare Rummel's formulation vs the null space approach.
% Test in small radius orbits using Linear Leas Squares.
% Author: Sergio Coll
% Date: 09/28/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % path = "HARMCOEFS_EROS_CD_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % n_max  = 4;
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
r      = 0.36E3;
% % r      = 24E3;
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
f = 1/10; dt = 1/f;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);
At = t(2) - t(1);

% position, orientation and angulat velocity errors
sigmaP = 0.01;                    % orbit uncertainty   [m]
sigmaE = 10 * pi / (3600 * 180);  % attitude errors     [rads]
sigmaW = 1E-3;                    % attitude errors     [deg/sqrt(hr)]

% state errors
sigmaW        = sigmaW / sqrt(At) * pi / (180*3600);                        % [rad/s]
deltaPos      = normrnd(0, sigmaP, [3, Nt]) + [1;1;1].*0.01;                % [m]
deltaAtt      = normrnd(0, sigmaE, [3, Nt]) + [1;6;10].*(pi/(180*3600));    % [rad]
deltaOmega    = normrnd(0, sigmaW, [3, Nt]);                                % [rad/s]

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(5+Nc+Ns,5+Nc+Ns), [(5+Nc+Ns)^2, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
dX = ones(length(X)-1, 1) * 1E-1;

% Gravity estimation
sigma  = 1E-11 * sqrt(f);   % [10 mE / sqrt(Hz)]
R  =  diag([1,1,1,1,1,1])*sigma^2;
STD = sqrt(diag(R)');
mu = zeros(length(STD), 1);
Nm = length(R);

% create measurement noise
noise0 = zeros(9, Nt);
noise = normrnd(repmat(mu, 1, Nt), repmat(STD(:), 1, Nt));

Iner = [580, 0, 0;...
    0,640,0;...
    0,0,1100];
Ax1_NSM = 0; Nx1_NSM = 0;
Ax2_NSM = 0; Nx2_NSM = 0;
Ax3_NSM = 0; Nx3_NSM = 0;

Ax1_PEP = 0; Nx1_PEP = 0;
Ax2_PEP = 0; Nx2_PEP = 0;
Ax3_PEP = 0; Nx3_PEP = 0;

ns  = 6;
Ax2_SRIF = 0; Nx2_SRIF = 0;
Ax3_SRIF = 0; Nx3_SRIF = 0;
CN_NSM = zeros(1, Nt); CN_PEP = zeros(1,Nt);
for j = 1:Nt
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    rn_ACI = rn(:, j);
    
    % Inertially fixed
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

    w = [sqrt(GM/vecnorm(rn_ACI)^3);0;0];

    % computed meas. & partials
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(5, 2:end);...
         Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
    
     % compute Point mass approx error
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM,...
        rn_ACI, ACAF_ACI, ACAF_BODY);
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = ...
        compute_angularDyadPartials_v2(w, Iner);
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];
    Hrot_grad = [Hrot_grad(1:3, :); Hrot_grad(5:6, :);Hrot_grad(9, :)];
    Hrot_omega_dyad = [Hrot_omega_dyad(1:3, :); ...
        Hrot_omega_dyad(5:6, :);Hrot_omega_dyad(9, :)];
    H_omegaDot_dyad = [H_omegaDot_dyad(1:3, :); ...
        H_omegaDot_dyad(5:6, :);H_omegaDot_dyad(9, :)];

    % nuisance parameter partials (3 cases)
    Hp_1 = Hpos;                                                % well conditioned
    Hp_2 = [Hpos, Hrot_grad];                                   % ill-conditioned
    Hp_3 = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];        % ill-conditioned

    % generate residuals (linear approximation. Assume constant errors)
    dY_1 = Hc * dX + Hp_1 * deltaPos(:, j)                    + noise(:, j);
    dY_2 = Hc * dX + Hp_2 * [deltaPos(:, j);deltaAtt(:,j)]    + noise(:, j);
    dY_3 = Hc * dX + Hp_3 * [deltaAtt(:, j);deltaOmega(:, j)] + noise(:, j);

    % compute null space
    C_1 = null(Hp_1');
    [~,~,D] = svd(Hp_2');
    C_2 = D(:, 5);
    [~,~,D] = svd(Hp_3');
    C_3 = D(:, 5:6);

    r_1 = C_1' * R * C_1;
    r_2 = C_2' * R * C_2;
    r_3 = C_3' * R * C_3;

    % compute projector
    PA1 = Hp_1 * pinv(Hp_1' * (R\Hp_1)) * (Hp_1'/(R));
    V_1 = eye(Nm,Nm) - PA1;
    PA2 = Hp_2 * pinv(Hp_2' * (R\Hp_2)) * (Hp_2'/(R));
    V_2 = eye(Nm,Nm) - PA2;
    PA3 = Hp_3 * pinv(Hp_3' * (R\Hp_3)) * (Hp_3'/(R));
    V_3 = eye(Nm,Nm) - PA3;

    % Solve normal eq (case 1)
    ax1_NSM = (C_1' * Hc)' * inv(r_1) * (C_1' * Hc);
    nx1_NSM = (C_1' * Hc)' * (r_1 \ (C_1' * dY_1));

    ax1_PEP = (V_1 * Hc)' * inv(R) * (V_1 * Hc);
    nx1_PEP = (V_1 * Hc)' * (R \ (V_1 * dY_1));

    Ax1_NSM = Ax1_NSM + ax1_NSM; Nx1_NSM = Nx1_NSM + nx1_NSM;
    Ax1_PEP = Ax1_PEP + ax1_PEP; Nx1_PEP = Nx1_PEP + nx1_PEP;
    
    CN_NSM(j) = cond(Ax1_NSM); CN_PEP(j) = cond(Ax1_PEP);

    % Solve normal eq (case 2)
    ax2_NSM = (C_2' * Hc)' * inv(r_2) * (C_2' * Hc);
    nx2_NSM = (C_2' * Hc)' * inv(r_2) * (C_2' * dY_2);

    ax2_PEP = (V_2 * Hc)' * (R \ (V_2 * Hc));
    nx2_PEP = (V_2 * Hc)' * (R \ (V_2 * dY_2));

% %     ax2_PEP = P2 * (Hc' * (R\Hc));
% %     nx2_PEP = P2 * (Hc' * (R\dY_2));

    Ax2_NSM = Ax2_NSM + ax2_NSM; Nx2_NSM = Nx2_NSM + nx2_NSM;
    Ax2_PEP = Ax2_PEP + ax2_PEP; Nx2_PEP = Nx2_PEP + nx2_PEP;

    % Solve normal eq (case 3)
    ax3_NSM = (C_3' * Hc)' * inv(r_3) * (C_3' * Hc);
    nx3_NSM = (C_3' * Hc)' * inv(r_3) * (C_3' * dY_3);

    ax3_PEP = (V_3 * Hc)' * (R \ (V_3 * Hc));
    nx3_PEP = (V_3 * Hc)' * (R \ (V_3 * dY_3));

    Ax3_NSM = Ax3_NSM + ax3_NSM; Nx3_NSM = Nx3_NSM + nx3_NSM;
    Ax3_PEP = Ax3_PEP + ax3_PEP; Nx3_PEP = Nx3_PEP + nx3_PEP;
end

% show condition number ratio NSM / PEP
disp('NSM vs PEP condition number ratio, case 1 = ' +...
    string(cond(Ax1_NSM)/cond(Ax1_PEP)));
disp('NSM vs PEP condition number ratio, case 2 = ' +...
    string(cond(Ax2_NSM)/cond(Ax2_PEP)));
disp('NSM vs PEP condition number ratio, case 3 = ' +...
    string(cond(Ax3_NSM)/cond(Ax3_PEP)));

Nx  = length(dX);

% plot condition number
figure()
semilogy(t, CN_NSM, t, CN_PEP, 'LineWidth', 2)
legend('NSM', 'PEP')
title('Condition number Case 1')
grid on;

% Solve Normal equations (case 1)
dX_NSM = Ax1_NSM \ Nx1_NSM;
dX_PEP = Ax1_PEP \ Nx1_PEP;

Xerr_NSM = abs(dX_NSM - dX);
Xerr_PEP = abs(dX_PEP - dX);

sigma_NSM = sqrt(diag(inv(Ax1_NSM)));
sigma_PEP = sqrt(diag(inv(Ax1_PEP)));
deltaErr_1 = abs(Xerr_NSM - Xerr_PEP)./abs(Xerr_PEP);

figure();
semilogy(1:Nx, 3.*sigma_NSM, 'LineWidth', 2, 'Color', 'b');
hold all;
semilogy(1:Nx, 3.*sigma_PEP, 'LineWidth', 2, 'Color', 'r');
semilogy(1:Nx, Xerr_NSM, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
semilogy(1:Nx, Xerr_PEP, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
title('Case 1')
grid on;

% Solve Normal equations (case 2)
dX_NSM = Ax2_NSM \ Nx2_NSM;
dX_PEP = Ax2_PEP \ Nx2_PEP;

Xerr_NSM = abs(dX_NSM - dX);
Xerr_PEP = abs(dX_PEP - dX);

sigma_NSM = sqrt(diag(inv(Ax2_NSM)));
sigma_PEP = sqrt(diag(inv(Ax2_PEP)));

deltaErr_2 = abs(Xerr_NSM - Xerr_PEP)./abs(Xerr_PEP);

figure();
semilogy(1:Nx, 3.*sigma_NSM, 'LineWidth', 2, 'Color', 'b');
hold all;
semilogy(1:Nx, 3.*sigma_PEP, 'LineWidth', 2, 'Color', 'r');
semilogy(1:Nx, Xerr_NSM, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
semilogy(1:Nx, Xerr_PEP, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
title('Case 2')
grid on;

% Solve Normal equations (case 3)
dX_NSM = Ax3_NSM \ Nx3_NSM;
dX_PEP = Ax3_PEP \ Nx3_PEP;

Xerr_NSM = abs(dX_NSM - dX);
Xerr_PEP = abs(dX_PEP - dX);

sigma_NSM = sqrt(diag(inv(Ax3_NSM)));
sigma_PEP = sqrt(diag(inv(Ax3_PEP)));
deltaErr_3 = abs(Xerr_NSM - Xerr_PEP)./abs(Xerr_PEP);

figure();
semilogy(1:Nx, 3.*sigma_NSM, 'LineWidth', 2, 'Color', 'b');
hold all;
semilogy(1:Nx, 3.*sigma_PEP, 'LineWidth', 2, 'Color', 'r');
semilogy(1:Nx, Xerr_NSM, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
semilogy(1:Nx, Xerr_PEP, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
title('Case 3')
grid on;

% plot error difference
[~, xvals] = orderValues(deltaErr_3, n_max);
figure();
semilogy(xvals, deltaErr_1.*100, 'LineWidth', 2, 'Color', 'b');
hold on;
semilogy(xvals, deltaErr_2.*100, 'LineWidth', 2, 'Color', 'r');
xlabel('Degree')
ylabel('[%]')
grid on;
title('Error ratio for the NSM and ANSM to Pre-elimination')
legend('Trajectory errors', 'Trajectory + atittude errors')

%% FUNCTIONS
function [yvals, xvals] = orderValues(vals, n_max)
    [Nc, ~, Ncs] = count_num_coeff(n_max); 
    C_vals = vals(1:Nc-1);
    S_vals = vals(Nc:Ncs-1);

    % get values starting at n = 2
    yvals = []; xvals = [];
    for n = 2:n_max
        [Nc, Ns, ~] = count_num_coeff(n); 

        C_maxIdx = Nc-1;
        C_minIdx = (Nc-1) - n;

        S_maxIdx = Ns;
        S_minIdx = Ns - (n-1); 

        vec = [C_vals(C_minIdx:C_maxIdx); S_vals(S_minIdx:S_maxIdx)];
        yvals = [yvals; vec];

        Nvals = length(vec);
        Nmin = n;
        Nmax = n + 1;
        xvals    =  [xvals, linspace(Nmin, Nmax-0.1, Nvals)];
    end
end