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

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
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
Nt = length(t);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(5+Nc+Ns,5+Nc+Ns), [(5+Nc+Ns)^2, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% gradiometer noise level
sigma  = 1;
R      = diag([1,1,1,1,1,1])*sigma^2;
deltaR = 1;
STD    = sqrt(diag(R)');
mu     = zeros(length(STD), 1);
Nm     = length(R);

% output variables
Ax0_NSM_pos  = 0; Ax0_NSM_att  = 0;
Ax1_NSM_pos  = zeros(Nc+Ns-1, Nc+Ns-1, length(R)); Ax1_NSM_att  = Ax1_NSM_pos;
for j = 1:Nt
    disp(j/Nt * 100)
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % Inertially-fixed
    rn_ACI = rn(:, j);
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

     % computed meas. & partials
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, ...
        [rn(:, j)', vn(:, j)'],zeros(9,1), Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
        Hc_BODY(5, 2:end); Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];

    % compute 1st order position partials
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
        rn_ACI, ACAF_ACI, ACAF_BODY);
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

    % compute attitude partials
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    Hrot_grad = [Hrot_grad(1:3, :); Hrot_grad(5:6, :);Hrot_grad(9, :)];
    
    % Information matrix NSM (position elimination)
    C = null(Hpos');
    H = C' * Hc; r = C' * R * C;
    Ax0_NSM_pos = Ax0_NSM_pos + (H'*(r\H));

    for k = 1:length(R)
        Rp = R;
        Rp(k,k) = R(k,k) + deltaR;
        r = C' * Rp * C;
        Ax1_NSM_pos(:, :, k) = Ax1_NSM_pos(:, :, k) + (H'*(r\H));
    end

    % Information matrix NSM (attitude elimination)
    C = null(Hrot_grad');
    H = C' * Hc; r = C' * R * C;
    Ax0_NSM_att = Ax0_NSM_att + (H'*(r\H));

    for k = 1:length(R)
        Rp = R;
        Rp(k,k) = R(k,k) + deltaR;
        r = C' * Rp * C;
        Ax1_NSM_att(:, :, k) = Ax1_NSM_att(:, :, k) + (H'*(r\H));
    end
end
ll= ["\delta \sigma^2_{xx} = 0.1", "\delta \sigma^2_{xy} = 0.1",...
    "\delta \sigma^2_{xz} = 0.1", "\delta \sigma^2_{yy} = 0.1",...
    "\delta \sigma^2_{yz} = 0.1", "\delta \sigma^2_{zz} = 0.1"];
figure()
hold all;
for k = 1:length(R)
    val  = squeeze(Ax1_NSM_pos(:, :, k));
    error = diag(Ax0_NSM_pos - val)./diag(Ax0_NSM_pos).*100;
    [RMS_error] = compute_RMS(error, n_max);
    plot(1:n_max, RMS_error, 'LineStyle','-', 'LineWidth', 2, ...
        'Marker','diamond', 'MarkerFaceColor','auto');
    hold on
end
xlabel('Degree')
ylabel('[%]')
xticks(1:n_max);          % positions
title('Sensitivity to noise mapping for ' + string(target_body))
legend(ll, 'NumColumns', 3);
grid on;

figure()
hold all;
for k = 1:length(R)
    val  = squeeze(Ax1_NSM_att(:, :, k));
    error = diag(Ax0_NSM_att - val)./diag(Ax0_NSM_att).*100;
    [RMS_error] = compute_RMS(error, n_max);
    plot(1:n_max, RMS_error, 'LineStyle','--', 'LineWidth', 2, ...
        'Marker','square', 'MarkerFaceColor','auto');
    hold on
end
xlabel('Coefficient Number')
ylabel('[%]')
xticks([0 2 3 4 5 6]);          % positions
title('Sensitivity to noise variance increament. Attitude NS')
legend(ll, 'NumColumns', 2);
grid on;