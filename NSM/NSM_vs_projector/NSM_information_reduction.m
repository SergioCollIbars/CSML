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

% traget body
target_body = "Bennu";  % options: Bennu / Eros

if(target_body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
    W = 4.06130329511851E-4;          % Rotation ang. vel   [rad/s]
    W0 = 0;                           % Initial asteroid longitude
    RA = deg2rad(86.6388);            % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);          % Declination         [rad]
    r = 0.36E3;                       % orbit radius        [m]
elseif(target_body == "Eros")
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    n_max  = 10;
    normalized = 1;
    GM =  459604.431484721;          % Point mass value    [m^3/s^2]
    W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
    W0 = 0;                          % Initial asteroid longitude
    RA = deg2rad(11.363);            % Right Ascension     [rad]
    DEC = deg2rad(17.232);           % Declination         [rad]
    r = 24E3;                        % orbit radius        [m]
end
poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

Cnm = Cnm(1:n_max+1, 1:n_max+1);
Snm = Snm(1:n_max+1, 1:n_max+1);

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
phi    = pi/2 - deg2rad(0);
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 3;
f   = 1/10;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(5+Nc+Ns,5+Nc+Ns), [(5+Nc+Ns)^2, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% gradiometer noise level
sigma  = 10 * sqrt(f); % [mE]
sigma  = 1;
R      = diag([1,1,1,1,1,1])*sigma^2;
STD    = sqrt(diag(R)');
mu     = zeros(length(STD), 1);
Nm     = length(R);

noise0 = zeros(9, Nt);
noise = normrnd(repmat(mu, 1, Nt), repmat(STD(:), 1, Nt));

% plot 3D trajectory
figure()
subplot(1, 2, 1)
plot3(rn(1, :), rn(2, :), rn(3, :), 'LineWidth', 2)
axis equal;
grid on;
subplot(1, 2, 2)
plot(t, vecnorm(rn) - Re)

% output variables
Ax_original = 0;
Ax_NSM_pos  = 0; Ax_NSM_att  = 0;

H_stack     = nan(6 * Nt, Ncs - 1);
G1_stack    = nan(6 * Nt, 3); G2_stack = G1_stack;
fprintf('  0%%');
for j = 1:Nt
    pct = 100 * j / Nt;
    fprintf('\b\b\b\b%3.0f%%', pct);

    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % Inertially-fixed
    rn_ACI    = rn(:, j);
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

     % computed meas. & partials
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, ...
        [rn(:, j)', vn(:, j)'],noise0, Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
        Hc_BODY(5, 2:end); Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)]./1E-12;

    % compute 1st order position partials
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
        rn_ACI, ACAF_ACI, ACAF_BODY);
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;

    % compute 2nd order position partials
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    Hrot_grad = [Hrot_grad(1:3, :); Hrot_grad(5:6, :);Hrot_grad(9, :)]./1E-12;
    
    % original information matrix
    Ax_original = Ax_original + (Hc'*(R\Hc));

    maxInd = 6 * j; minInd = maxInd - 5;
    H_stack(minInd:maxInd, :) = Hc;
    G1_stack(minInd:maxInd, :) = Hpos;
    G2_stack(minInd:maxInd, :) = Hrot_grad;

    % Information matrix NSM (position elimination)
    C = null(Hpos');
    H = C' * Hc; r = C' * R * C;
    Ax_NSM_pos = Ax_NSM_pos + (H'*(r\H));

    % Information matrix NSM (attitude elimination)
    C = null(Hrot_grad');
    H = C' * Hc; r = C' * R * C;
    Ax_NSM_att = Ax_NSM_att + (H'*(r\H));
end
% % P_G             = eye(Nt*6) - G1_stack * (pinv(G1_stack' * G1_stack) * G1_stack');
% % Ax_original_pos = H_stack' * P_G * H_stack;

Ht              = [H_stack, G1_stack];
F               = (sigma^2).*inv(Ht' * Ht);
Ax_original_pos = inv(F(1:Ncs-1, 1:Ncs-1));


% % P_G             = eye(Nt*6) - G2_stack * ((G2_stack' * G2_stack) \ G2_stack');
% % Ax_original_att = H_stack' * P_G * H_stack;

Ht              = [H_stack, G2_stack];
F               = (sigma^2).*inv(Ht' * Ht);
Ax_original_att = inv(F(1:Ncs-1, 1:Ncs-1));

% Error fraction for position
error_pos_std   = sqrt(diag(inv(Ax_NSM_pos))./diag(inv(Ax_original_pos)));
error_pos_info  = (diag(Ax_NSM_pos)./diag(Ax_original_pos));

% Error fraction for attitude
error_att_std   = sqrt(diag(inv(Ax_NSM_att))./diag(inv(Ax_original_att)));
error_att_info  = (diag(Ax_NSM_att)./diag(Ax_original_att));

% correlations
P_NSM_att       = inv(Ax_NSM_att);
P_NSM_pos       = inv(Ax_NSM_pos);
P_orig_pos      = inv(Ax_original_pos);
P_orig_att      = inv(Ax_original_att);

F0 = abs(triu(corr(P_orig_pos), 1));
F1 = abs(triu(corr(P_orig_att), 1));
F2 = abs(triu(corr(P_NSM_pos), 1));
F3 = abs(triu(corr(P_NSM_att), 1));

figure()
subplot(1, 2, 1)
heatmap(F0, 'Colormap',turbo)
title('Original Pos. Cov.')
subplot(1, 2, 2)
heatmap(F2, 'Colormap',turbo)
title('NSM position Cov.')
figure()
subplot(1, 2, 1)
heatmap(F1, 'Colormap',turbo)
title(' Original Att. Cov.')
subplot(1, 2, 2)
heatmap(F3, 'Colormap',turbo)
title('NSM attitude Cov.')

% RMS values STD ratio
[error_pos_RMS]        = compute_RMS(error_pos_std, n_max);
[yvals_pos, xvals_pos] = orderValues(error_pos_std, n_max);
[error_att_RMS]        = compute_RMS(error_att_std, n_max);
[yvals_att, xvals_att] = orderValues(error_att_std, n_max);

tt = 'SH STD increment for ' + string(target_body);
yy = '\Delta \sigma_{C,S}';
plot_results(n_max, error_pos_RMS, error_att_RMS, ...
    yvals_pos, yvals_att, xvals_pos, xvals_att, tt, yy)

% Info ratio
[error_pos_RMS]        = compute_RMS(error_pos_info, n_max);
[yvals_pos, xvals_pos] = orderValues(error_pos_info, n_max);
[error_att_RMS]        = compute_RMS(error_att_info, n_max);
[yvals_att, xvals_att] = orderValues(error_att_info, n_max);

tt = 'Information retained for ' + string(target_body);
yy = '[%]';
plot_results(n_max, error_pos_RMS.*100, error_att_RMS.*100, ...
    yvals_pos.*100, yvals_att.*100, xvals_pos, xvals_att, tt, yy);


%% FUNCTIONS
function [] = plot_results(n_max, error_pos_RMS, error_att_RMS, ...
    yvals_pos, yvals_att, xvals_pos, xvals_att, tt, yy)
    figure()
    plot((1:n_max), error_pos_RMS, 'LineWidth', 2, ...
        'Marker','square', 'Color', 'b')
    hold all;
    plot((1:n_max), error_att_RMS, 'LineWidth', 2, 'Marker', 'square', ...
        'Color','r')
    plot(xvals_pos, yvals_pos, 'LineStyle', 'none', 'Marker','*', ...
        'Color','b');
    plot(xvals_att, yvals_att, 'LineStyle', 'none', 'Marker','*', ...
        'Color','r');
    xlabel('Degree')
    ylabel(yy)
    title(tt)
    legend('NSM position', 'NSM attitude');
    grid on; xticks(1:n_max);
end