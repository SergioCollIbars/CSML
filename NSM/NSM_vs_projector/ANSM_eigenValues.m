clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%%              ANSM PLOT EIGENVALUES
% Description: Plot the smallest eigenvalue for the ANSM
% Author: Sergio Coll
% Date: 10/27/24

target_body = "Eros";  % options: Bennu or Eros

if(target_body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
    W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
    W0 = 0;                   % Initial asteroid longitude
    RA = deg2rad(86.6388);    % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);  % Declination         [rad]
    r = 0.36E3;               % orbit radius        [m]
elseif(target_body == "Eros")
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    n_max  = 4;
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
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 3;
f   = 1/10;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);
At = t(2) - t(1);

% Integrate trajectory (true trajectory)
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 14, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% perturb nominal coefficient
[X]   = mat2list(Cnm, Snm, Nc, Ns);
Iner = [580, 0, 0;...
    0,640,0;...
    0,0,1100];

% compute eigenvalues
eigenVal = ones(6, Nt) * NaN;
delta_n  = ones(3, 1);
for j = 1:Nt
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    rn_ACI = rn(:, j);
    
    % Inertially fixed
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

    w = [0;0;sqrt(GM/vecnorm(rn_ACI)^3)];

    % computed meas. & partials
    [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
        [rn(:, j)', vn(:, j)'], zeros(9,1), Cnm, Snm);
    Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
    Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end);...
        Hc_BODY(5, 2:end); Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
    
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
    scale = 1E-9;                                                      % to Eotvos
    Hn = [Hrot_grad, (Hrot_omega_dyad+H_omegaDot_dyad).*At]./(scale);  % (ill-conditioned)
% %     Hn = Hrot_grad./scale;

    % SVD decomposition
    [S, V, D] = svd(Hn');
    eigenVal(:, j) = diag(V);
% %     F = sum(V)';
% % 
% %     % bound
% %     a = abs(D' * Hn * delta_n);
% %     b = F.*vecnorm(delta_n);
end

figure()
tt = ["\lambda_1","\lambda_2", "\lambda_3",...
    "\lambda_4","\lambda_5","\lambda_6" ];
semilogy(t./T, eigenVal, 'LineWidth', 2)
grid on;
xlabel('Orbit rev.');
ylabel('Eotvos / rad');
xticks(1:T);
legend(tt, 'NumColumns', 2);
title('Nuisance eigenvalues for ' + string(target_body))
