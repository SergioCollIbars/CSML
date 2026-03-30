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
saveData    = 0;
loadData    = 0;
mc_num      = 1;

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
    n_max  = 4;
    normalized = 1;
    GM =  459604.431484721;          % Point mass value    [m^3/s^2]
    W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
    W0 = 0;                          % Initial asteroid longitude
    RA = deg2rad(11.363);            % Right Ascension     [rad]
    DEC = deg2rad(17.232);           % Declination         [rad]
    r = 36E3;                        % orbit radius        [m]
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
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% scale states
scale = 1;
r0 = r0 * scale;         % [-]
v0 = v0 * scale;         % [- / s]
GM = GM * (scale)^3;     % [-^3/s^2]
Re = Re * scale;         % [-]
asterParams(1) = GM; asterParams(2) = Re;

if(loadData == 1)
    load("Initial_conditions.mat");
    r0 = X_true_final(1:3)'; v0 = X_true_final(4:6)';
    t = t + t_final;
end

% Integrate trajectory (true trajectory)
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 14, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);

% perturb nominal coefficient
[X]   = mat2list(Cnm, Snm, Nc, Ns);
sigmaCoeff = 1E-3;

% Integrate trajectory (nominal trajectory)
if(loadData == 1)
    X0 = X_estim_final(end-5:end);
    [Cp, Sp] = list2mat(n_max, Nc, Ns, X_estim_final(1:Nc+Ns));
    P0_LS = P_final;
    P0_NSM = P_final(1:Ncs-1, 1:Ncs-1);
end

% Gradiometer measurement noise
sigma  = 1E-12 * sqrt(f);
R  =  diag([1,1,1,1,1,1])*sigma^2;
RR =  diag([1,1,1,1,1,1,1,1,1])*sigma^2;
STD = sqrt(diag(RR)');
mu = zeros(length(STD), 1);
Nm = length(RR);

noise0 = zeros(9, Nt);

dX1 = ones(Ncs, mc_num) * NaN; dX3 = ones(Ncs+6, mc_num) * NaN;
for Mc = 1:mc_num
    disp('Mc run = ' + string(Mc))
    % create gradiometer measurements
    Y_true = ones(9, length(t));
    noise = normrnd(repmat(mu, 1, Nt), repmat(STD(:), 1, Nt));
    for k = 1:length(t)
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        r = state_t(k, 1:3)'; v = state_t(k, 4:6)';
        [Y_ACI, ~] = gradiometer_meas(t(k) ,[GM, Re, 14, normalized],...
            ACAF_ACI, [r', v'], zeros(9, 1), Cnm, Snm);
    
        Y_true(:, k) = Y_ACI + noise(:, k);
    end

    % create initial conditions
    deltaR = normrnd(0, sigmaP, [3, 1]);
    deltaX = normrnd(0, sigmaCoeff, [Ncs, 1]);

    X0 = [r0+deltaR;v0];
    deltaX(1) = 0;
    [Cp, Sp] = list2mat(n_max, Nc, Ns, X+deltaX);
    P0_LS = diag([100*deltaX(2:end); 100.*ones(6, 1)]);
    P0_NSM = diag(100*deltaX(2:end));
    
    % generate initial nominal orbit
    STM0 = reshape(eye(6+Nc+Ns-1,6+Nc+Ns-1), [(6+Nc+Ns-1)^2, 1]);
    [~, state_n] = ode113(@(t, x) EoM(t, x, Cp, Sp, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 1), t, [X0;STM0], options);
    rn = state_n(:, 1:3)';
    vn = state_n(:, 4:6)';
    dr = state_t(:, 1:3)' - rn;
    
    [X_LS, P_LS] = solver_LS(t, X0(1:3), X0(4:6), Cp, Sp, asterParams, ...
        poleParams, P0_LS, Y_true, R);
    
    [X_NSM, P_NSM] = solver_NSM(t, X0(1:3), X0(4:6), Cp, Sp, asterParams, ...
        poleParams, P0_NSM, Y_true, R);
    
    P1  = P_NSM; P3  = P_LS;
    dX1(:, Mc) = X_NSM; dX3(:, Mc) = X_LS;
end

% plot condition number over time
figure()
semilogy(t./T, cond_CEP, t./T, cond_NSM, 'LineWidth', 2);
legend('Co-estimation: position, velocity + grav. field', 'NSM')
title('Condition number over time');
xlabel('Orbit revolution');
grid on;

% Compute correlation matrix
std_devs = sqrt(diag(P1));
R1 = P1./ (std_devs * std_devs');

std_devs = sqrt(diag(P3));
R3 = P3./ (std_devs * std_devs');

% Plot correlations
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

% Compute actual errors
Nx  = length(deltaX);

Xerr1 = abs(dX1 - X);
Xerr3 = abs(dX3(1:Nx, :) - X);
Perr  = abs(dX3(end-5:end, :) - [r0;v0]);

sigma1 = sqrt(diag(P1));
sigma3 = sqrt(diag(P3));

figure()
semilogy(1:3, 3.*sigma3(end-5:end-3), 'Color', 'r', 'LineWidth', 2);
hold all
for Mc = 1:mc_num
    semilogy(1:3, Perr(1:3, Mc), 'Color', 'r', 'Marker','*', ...
        'LineStyle','none');
end

figure()
semilogy(1:3, 3.*sigma3(end-2:end), 'Color', 'r', 'LineWidth', 2);
hold all
for Mc = 1:mc_num
    semilogy(1:3, Perr(4:6, Mc), 'Color', 'r', 'Marker','*', ...
        'LineStyle','none');
end

figure()
semilogy(1:Nx-1, 3.*sigma1, 'Color', 'b', 'LineWidth', 2);
hold all;
semilogy(1:Nx-1, 3.*sigma3(1:Nx-1), 'Color', 'r', 'LineWidth', 2);
for Mc = 1:mc_num
    semilogy(1:Nx-1, Xerr1(2:Nx, Mc), 'Color', 'b', 'Marker','*', ...
        'LineStyle','none')
    hold on;
    semilogy(1:Nx-1, Xerr3(2:Nx, Mc), 'Color', 'r', 'Marker','*', ...
        'LineStyle','none');
end
legend('\sigma NSM', '\sigma Co-estim', 'err NSM', 'err Co-estim')
grid on;

% plot error along trajectory
% % figure()
% % plot(t./T, Arv_t(1:3, :)./scale, 'LineWidth', 2);
% % xlabel('Orbit revolution')
% % ylabel('[m]')
% % title('Position error over time')
% % grid on;

%% FUNCTIONS

function [X_final, P_final] = solver_LS(t, rn0, vn0, Cnm, Snm, asterParams, ...
    poleParams, P0, Y_true, R)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    Nx = Ncs + 5;

    % integrator options
    STM0 = reshape(eye(Nx), [Nx*Nx, 1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    
    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns);rn0;vn0];

    count = 0;
    iterMax = 5;
    err = 0;
    xnot_L = zeros(Ncs-1+6, 1);
    while (count < iterMax) && (err < 1E3)
        % integrate state
        [~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
            W0, W, RA, DEC, 1), t, [rn0;vn0;STM0], options);
        rn = state_n(:, 1:3)';
        vn = state_n(:, 4:6)';

        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
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
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
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
            H = [Hc, Hpos, zeros(6, 3)] * PHI;

            % prefit residuals
            dY = Y_true(:, j) - Y_ACI;
            dy = [dY(1:3);dY(5:6);dY(9)];

            % solve normal equations
            ax = H' * (R\H);
            nx = H' * (R\dy);
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
        rn0 = X_new(Ncs+1:Ncs+3); vn0 = X_new(Ncs+4:end);
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  LS update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    P_final = inv(Ax_L);
end


function [X_final, P_final] = solver_NSM(t, rn0, vn0, Cnm, Snm, asterParams, ...
    poleParams, P0, Y_true, R)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    Nx = Ncs + 5;

    % integrator options
    STM0 = reshape(eye(Nx), [Nx*Nx, 1]);
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    
    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 5;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    while (count < iterMax) && (err < 1E3)
        [~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
            W0, W, RA, DEC, 1), t, [rn0;vn0;STM0], options);
        rn = state_n(:, 1:3)';
        vn = state_n(:, 4:6)';

        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
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

            % prefit residuals
            dY = Y_true(:, j) - Y_ACI;
            dy = C' * [dY(1:3);dY(5:6);dY(9)];

            % solve normal equations
            ax = H' * (r\H);
            nx = H' * (r\dy);
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  NSM update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    P_final = inv(Ax_L);
end