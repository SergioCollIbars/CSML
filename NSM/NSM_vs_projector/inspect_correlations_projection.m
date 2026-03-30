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
r      = 0.35E3;
% % r      = 24E3;
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

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

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
sigma  = 1E-12;
% % R = diag([1,1,1,1,1,1,1,1,1])*sigma^2;
% % R  = diag([100.7,1,6,7.6,2,4000])*sigma^2;
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
Ax1 = 0; Nx1 = 0;
Ax2 = 0; Nx2 = 0;
Ax3 = 0; Nx3 = 0;
cond_CEP = zeros(1, Nt); prefit1 = zeros(3, Nt);
cond_NSM = zeros(1, Nt); prefit2 = zeros(6, Nt);
H1_store = zeros(3*Nt, Ncs-1); 
H2_store = zeros(6*Nt, Ncs-1); 
for j = 1:Nt
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % RTN rotation matrix
    ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    rn_ACI = rn(:, j);
    
% %     ACAF_BODY  = ACAF_ACI * ACI_RTN;
% %     B_ACI = ACI_RTN';

    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

    w = [0;0;sqrt(GM/vecnorm(rn_ACI)^3)];

    % computed meas.
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
% %     Hc = Hc_BODY;
% % 
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); Hc_BODY(5, 2:end);...
         Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
% %      Hc = [Hc_BODY(1, :); Hc_BODY(4, :); Hc_BODY(7, :); Hc_BODY(5, :);...
% %          Hc_BODY(8, :); Hc_BODY(9, :)];
% %     Hc = [Hc_BODY(1, :); Hc_BODY(4, :); Hc_BODY(7, :); Hc_BODY(2, :); Hc_BODY(5, :);...
% %                   Hc_BODY(8, :);  Hc_BODY(3, :); Hc_BODY(6, :); Hc_BODY(9, :)];
    
     % compute Point mass approx error
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, ACAF_BODY);
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = compute_angularDyadPartials_v2(w, Iner);
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];
    Hrot_grad = [Hrot_grad(1:3, :); Hrot_grad(5:6, :);Hrot_grad(9, :)];
    Hrot_omega_dyad = [Hrot_omega_dyad(1:3, :); Hrot_omega_dyad(5:6, :);Hrot_omega_dyad(9, :)];
    H_omegaDot_dyad = [H_omegaDot_dyad(1:3, :); H_omegaDot_dyad(5:6, :);H_omegaDot_dyad(9, :)];

    % % Hrot = [Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];
    % % Hrot = [Hpos, Hrot_grad, Hrot_omega_dyad+H_omegaDot_dyad];
    % % Hrot = [Hpos, Hrot_grad];
     Hrot = Hpos;

    % generate residuals
    dY = Hc * dX + Hrot * Ar + noise(:, j);
% %         dY = Hc * dX + Hrot * Ar + noise(:, j);
    PHI = reshape(state_t(j, 7:end), [5+Nc+Ns, 5+Nc+Ns]);
% %     F = PHI * [dX; Ar; zeros(3, 1)];
% %     Ari = [F(Nc+Ns:Nc+Ns+2)];
% %     dY = Hc * dX + Hrot * Ari + noise(:, j);
% %     dY = Hc * dX;

    % compute null space & compute error
    C = null(Hrot');
% %     [S,v,D] = svd(Hrot');
% %     C = D(:, 5:6);
    r = C' * R * C;
    err1 =  (C'* (dY - Hc * dX))./((C' * Hc) * dX);
    prefit1(:, j) = C' * dY;
    maxPos = 3 * j; minPos = maxPos- 2;
    H1_store(minPos:maxPos, :) = C' * Hc;

    % compute projector
    PA1 = Hrot * inv(Hrot' * (R\Hrot)) *(Hrot'/(R));
    V = eye(Nm,Nm) - PA1;

    err2 = (V * (dY - Hc * dX))./(V * Hc * dX);
    prefit2(:, j) = V * dY;
    maxPos = 6 * j; minPos = maxPos- 5;
    H2_store(minPos:maxPos, :) = V * Hc;

    % investigate correlations
    H = [Hc, Hrot, zeros(Nm, 3)] * PHI;
    % % B = [Hc, Hpos, zeros(6, 3)] * PHI(2:end, 2:end);
    % % H = [B, Hrot_grad];
    % % H = [Hc, Hpos, zeros(Nm, 3), Hrot_grad] * PHI(2:end, 2:end);
    % % H = [Hc, Hpos, zeros(Nm, 3)] * PHI;
    % % H = [Hc, Hrot];
    % % H = [Hc, Hrot; zeros(1, 46), B, zeros(1, 6)];
    % % H = Hc;

    % Solve normal eq
    ax1 = (C' * Hc)' * inv(r) * (C' * Hc);
    nx1 = (C' * Hc)' * inv(r) * (C' * dY);

    ax2 = (V * Hc)' * inv(R) * (V * Hc);
    nx2 = (V * Hc)' * inv(R) * (V * dY);

    ax3 = (H)' * inv(R) * (H);
    nx3 = (H)' * inv(R) * (dY);

    Ax1 = Ax1 + ax1; Nx1 = Nx1 + nx1;
    Ax2 = Ax2 + ax2; Nx2 = Nx2 + nx2;
    Ax3 = Ax3 + ax3; Nx3 = Nx3 + nx3;

    cond_NSM(j) = cond(Ax1);
    cond_CEP(j) = cond(Ax3);
end

P1 = inv(Ax1); 
P2 = inv(Ax2);
P3 = inv(Ax3);

% plot condition number over time
figure()
semilogy(t./T, cond_CEP, t./T, cond_NSM, 'LineWidth', 2);
legend('Co-estimation', 'NSM')
title('Condition number over time');
xlabel('Orbit revolution');
grid on;

% Compute correlation matrix
std_devs = sqrt(diag(P1));
R1 = P1./ (std_devs * std_devs');

std_devs = sqrt(diag(P2));
R2 = P2./ (std_devs * std_devs');

std_devs = sqrt(diag(P3));
R3 = P3./ (std_devs * std_devs');

% show condition number
disp('Condition number Ax1 = ' + string(cond(Ax1)))
disp('Condition number Ax2 = ' + string(cond(Ax2)))
disp('Ax1 vs Ax2 condition number ratio = ' + string(cond(Ax1)/cond(Ax2)));

% Plot it
lg = ["C_{20}", "C_{21}", "C_{22}", "C_{30}", "C_{31}", "C_{32}", "C_{33}", ...
    "S_{21}", "S_{22}", "S_{31}", "S_{32}", "S_{33}", "X", "Y", "Z"];
figure;
imagesc(abs(R3));
colormap jet;
colorbar;
axis equal tight;
xticks(1:size(R3,1));
yticks(1:size(R3,1));
xlabel('Variable Index');
ylabel('Variable Index');
title('Correlation Matrix. Co-Estimation method');

figure;
imagesc(abs(R1));
colormap jet;
colorbar;
axis equal tight;
xticks(1:size(R1,1));
yticks(1:size(R1,1));
xlabel('Variable Index');
ylabel('Variable Index');
title('Correlation Matrix. NSM');

figure;
imagesc(abs(R2));
colormap jet;
colorbar;
axis equal tight;
xticks(1:size(R2,1));
yticks(1:size(R2,1));
xlabel('Variable Index');
ylabel('Variable Index');
title('Correlation Matrix. Pre-Elimination of paramers');

Nx  = length(dX);

dX1 = P1 * Nx1;
dX2 = P2 * Nx2;
dX3 = P3 * Nx3;

Xerr1 = abs(dX1 - dX);
Xerr2 = abs(dX2 - dX);
Xerr3 = abs(dX3(1:Nx) - dX);

sigma1 = sqrt(diag(P1));
sigma2 = sqrt(diag(P2));
sigma3 = sqrt(diag(P3));

figure()
semilogy(1:Nx, 3.*sigma2, 1:Nx, 3.*sigma1, 1:Nx, 3.*sigma3(1:Nx), 'LineWidth', 2)
hold on;
semilogy(1:Nx, Xerr2, 1:Nx, Xerr1, 1:Nx, Xerr3, 'LineWidth', 2, 'LineStyle', '--')
legend('\sigma pre-elim', '\sigma NSM', '\sigma Co-estim', 'err pre-elim', 'err NSM', 'err Co-estim')
grid on;

% postfit residuals
residuals1 = prefit1.*0;
residuals2 = prefit2.*0;
for j = 1:Nt
    maxPos = 3*j; minPos = maxPos - 2;
    residuals1(:, j) = (prefit1(:, j) - H1_store(minPos:maxPos, :) * dX1);
        
    maxPos = 6*j; minPos = maxPos - 5;
    residuals2(:, j) =  (prefit2(:, j) - H2_store(minPos:maxPos, :) * dX2);
end
