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

target_body = "Bennu";  % options: Bennu or Eros

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
    r = 0.36E3;                  % orbit radius        [m]
    sigmaP = 0.03;
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
rev = 6;
f = 1/30;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% scale states
scale = 1;
r0 = r0 * scale;         % [-]
v0 = v0 * scale;         % [- / s]
GM = GM * (scale)^3;     % [-^3/s^2]
Re = Re * scale;         % [-]
asterParams(1) = GM; asterParams(2) = Re;

% Integrate trajectory (true trajectory)
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 14, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);

% perturb nominal coefficient
[X]   = mat2list(Cnm, Snm, Nc, Ns);
sigmaCoeff = 1E-3;

% Gradiometer measurement noise
sigma  = 10E-12 * sqrt(f);
R  =  diag([1,1,1,1,1,1])*sigma^2;
RR =  diag([1,1,1,1,1,1,1,1,1])*sigma^2;
STD = sqrt(diag(RR)');
mu = zeros(length(STD), 1);
Nm = length(RR);

noise0 = zeros(9, Nt);

X0       = [r0;v0]; scale = 1;
sigmaP   = .1 * scale;   % [m]
sigmaV   = 0.01 * scale; % [m/s]
deltaX   = normrnd(0, sigmaCoeff, [Ncs, 1]);
[Cp, Sp] = list2mat(n_max, Nc, Ns, X);
P0_LS_CS = diag([abs(deltaX(2:end));sigmaP.*ones(3, 1);sigmaV.*ones(3, 1)]);
P0_NSM   = diag(abs(deltaX(2:end)));
P0_LS    = diag(abs(deltaX(2:end)));

% generate initial nominal orbit
STM0 = reshape(eye(6+Nc+Ns-1,6+Nc+Ns-1), [(6+Nc+Ns-1)^2, 1]);
[~, state_n] = ode113(@(t, x) EoM(t, x, Cp, Sp, n_max, GM, Re, normalized, ...
W0, W, RA, DEC, 1), t, [X0;STM0], options);
rn = state_n(:, 1:3)';
vn = state_n(:, 4:6)';

[cond_CS] = solver_LS_CS(t, X0(1:3), X0(4:6), Cp, Sp, asterParams, ...
    poleParams, P0_LS_CS, R);

[cond_LS] = solver_LS(t, X0(1:3), X0(4:6), Cp, Sp, asterParams, ...
    poleParams, P0_LS, R);

[cond_NSM] = solver_NSM(t, X0(1:3), X0(4:6), Cp, Sp, asterParams, ...
    poleParams, P0_NSM, R);

% plot condition number over time
figure()
semilogy(t./T, cond_CS, 'Color' , 'b', 'LineWidth', 2); hold all;
semilogy(t./T, cond_NSM, 'Color' ,'g', 'LineWidth', 2);
semilogy(t./T, cond_LS, 'Color', 'magenta', 'LineWidth', 2);
legend('Gravity coefficients + position & velocity', ...
    'Gravity coefficients using NSM',...
    'Gravity coefficients (ideal model)')
title('Condition number over time');
xlabel('Orbit revolution');
grid on;

%% FUNCTIONS

function [cond_LS_CS] = solver_LS_CS(t, rn0, vn0, Cnm, Snm, asterParams, ...
    poleParams, P0, R)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [~, ~, Ncs] = count_num_coeff(n_max); 
    Nx = Ncs + 5;

    % integrator options
    STM0 = reshape(eye(Nx), [Nx*Nx, 1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    
    % white meas.
    Lk = chol(R, 'lower');       % Rk = Lk*Lk'

    % integrate state
    [~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 1), t, [rn0;vn0;STM0], options);
    rn = state_n(:, 1:3)';
    vn = state_n(:, 4:6)';

    Ax_L = inv(P0); cond_LS_CS = nan(1, Nt);
    for j = 1:Nt        
        % STM at current time
        PHI = reshape(state_n(j, 7:end), [Nx,Nx]);
        rn_ACI = rn(:, j);

        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        B_ACI     = eye(3,3);
        ACAF_BODY = ACAF_ACI * B_ACI';

        % computed meas & select measurements
        [~, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
            [rn(:, j)', vn(:, j)'], ...
                zeros(9,1), Cnm, Snm);
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
            Hc_BODY(5, 2:end);...
            Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
    
         % compute position and attitude partials
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
            rn_ACI, ACAF_ACI, ACAF_BODY);
        Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

        % visibility matrix
        H0 = [Hc, Hpos, zeros(6, 3)] * PHI;
        H  = [H0(:, 46:end),H0(:, 1:45)];
        
        R0 = chol(Ax_L);

        % SRIF factor via QR
        [~, RF] = qr([R0;Lk \ H], 0);      % A = Q*R, R upper triangular
        
        R = RF;

        % solve normal equations
% %         ax = H' * (R\H);

        ax = RF' * RF;
        
        Ax_L  =  ax;
        cond_LS_CS(j) = cond(R' * R);
    end
end

function [cond_LS] = solver_LS(t, rn0, vn0, Cnm, Snm, asterParams, ...
    poleParams, P0, R)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [~, ~, Ncs] = count_num_coeff(n_max); 
    Nx = Ncs + 5;

    % integrator options
    STM0 = reshape(eye(Nx), [Nx*Nx, 1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    

    % integrate state
    [~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 1), t, [rn0;vn0;STM0], options);
    rn = state_n(:, 1:3)';
    vn = state_n(:, 4:6)';

    Ax_L = inv(P0); cond_LS = nan(1, Nt);
    for j = 1:Nt        
        % STM at current time
        PHI = reshape(state_n(j, 7:end), [Nx,Nx]);
        rn_ACI = rn(:, j);

        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        B_ACI     = eye(3,3);
        ACAF_BODY = ACAF_ACI * B_ACI';

        % computed meas & select measurements
        [~, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
            [rn(:, j)', vn(:, j)'], ...
                zeros(9,1), Cnm, Snm);
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
            Hc_BODY(5, 2:end);...
            Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];

        % visibility matrix
        H = Hc;

        % solve normal equations
        ax = H' * (R\H);
        
        Ax_L  = Ax_L + ax;
        cond_LS(j) = cond(Ax_L);
    end
end

function [cond_NSM] = solver_NSM(t, rn0, vn0, Cnm, Snm, asterParams, ...
    poleParams, P0, R)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [~, ~, Ncs] = count_num_coeff(n_max); 
    Nx = Ncs + 5;

    % integrator options
    STM0 = reshape(eye(Nx), [Nx*Nx, 1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);

    [~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 1), t, [rn0;vn0;STM0], options);
    rn = state_n(:, 1:3)';
    vn = state_n(:, 4:6)';

    % integrate state
    Ax_L = inv(P0); cond_NSM = nan(1, Nt);
    for j = 1:Nt        
        % current position
        rn_ACI = rn(:, j);

        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        B_ACI     = eye(3,3);
        ACAF_BODY = ACAF_ACI * B_ACI';

        % computed meas & select measurements
        [~, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
            [rn(:, j)', vn(:, j)'], ...
                zeros(9,1), Cnm, Snm);
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
        Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
            Hc_BODY(5, 2:end);...
            Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
    
         % compute position and attitude partials
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
            rn_ACI, ACAF_ACI, ACAF_BODY);
        Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

        % Compute NS
        C = null(Hpos');

        % visibility matrix
        H = C' * Hc;
        r = C' * R * C;

        % solve normal equations
        ax = H' * (r\H);
        
        Ax_L  = Ax_L + ax;
        cond_NSM(j) = cond(Ax_L);
    end
end