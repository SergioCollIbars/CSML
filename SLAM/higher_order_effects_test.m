clear;
clc;
close all;

addpath("functions/"); addpath("data/");

set(0,'defaultAxesFontSize',16);

%%                  SIMULATION FRAMEWORK
% load meta data
metaData_path = "Metadata.txt";
mtd = readParams("data/"+metaData_path);

% start Simulation
[planetParams, poleParams, Cnm_t, Snm_t, ...
    initCond, t_range]                 = loadUniverse(metaData_path);
[instrumentParams]                     = loadInstrument(metaData_path);
[Nc, Ns, Ncs]  = count_num_coeff(planetParams(3)); 
X_true  = [mat2list(Cnm_t, Snm_t, Nc, Ns)];

[~, ~,  ~, ~, Qb, ~, ~, ~, ~] = ...
    loadFilterParams(metaData_path, planetParams, instrumentParams,...
    Cnm_t, Snm_t);

% integrate Trajectory
disp('Simulating trajectory ...')
tmin = t_range(1);
tmax = t_range(2);
t = linspace(tmin, tmax, tmax * instrumentParams(1, 5));

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
Nx      = 12 + Ncs; 
PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
[time, state] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
    Cnm_t, Snm_t), t, [initCond(1:3);initCond(4:6);...
    instrumentParams(:, 3);PHI0], options);

plot_trajectory(time, state)
disp('  DONE ...')

% trajectory error [m]
sigmaP    = 10;    Nt = length(time);                    
deltaR    = normrnd(0, sigmaP, [3, Nt]);

% compute S/C orientation
disp('Computing S/C orientation ...')
[BN_mat] = compute_orientation(time, state, instrumentParams);
disp('  DONE ...')

% generate Measurements
disp('Simulating measurements ...')
state_t = state; state_t(:, 1:3) = state_t(:, 1:3) + deltaR';
[Y_true, state_t] = compute_measurements(instrumentParams, planetParams, ...
    poleParams, time, state_t, Cnm_t, Snm_t, BN_mat);

plot_measurements(time, Y_true)
disp('  DONE ...')


%% ASSES HIGHER ORDER BOUNDS

% extract asteroid parameteres
GM = planetParams(1); Re = planetParams(2);
normalized = planetParams(4); 
t = time; Nt = length(time);
n_max = planetParams(3);

% extract asteroid pole parameters
W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
DEC = poleParams(4);

% nominal trajectory
rn = state(:, 1:3)'; vn = state(:, 4:6)';

% measurement noise
sigma_W = instrumentParams(1, 2) * sqrt(instrumentParams(1, 5));

Pc     = eye(3) * sigmaP^2;
s_mean = [Pc(1,1);Pc(1,2);Pc(1,3);Pc(2,2);Pc(2,3);Pc(3,3)];

% numerical error
fxx = deltaR(1, :).*deltaR(1, :);

% compute higher order effects
dY_HOT = nan(6, Nt); dY_FOT = nan(6, Nt); dY_SOT = nan(6, Nt);
mean_analy_HOT = nan(6, Nt);
bound  = nan(6, Nt);
for j = 1:Nt
    % current position
    rn_ACI = rn(:, j);
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI  =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    B_ACI     = eye(3,3);
    ACAF_BODY = ACAF_ACI * B_ACI';

    % computed measurements
    [Y_ACI, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI,...
        [rn(:, j)', vn(:, j)'], ...
            zeros(9,1), Cnm_t, Snm_t);
    y_nom = [Y_ACI(1:3); Y_ACI(5:6);Y_ACI(9)]./1E-12;

     % compute position and attitude partials
    [Hpos_tot] = compute_posPartials(n_max, normalized, Cnm_t, Snm_t, Re, GM, ...
        rn_ACI, ACAF_ACI, ACAF_BODY);
     Hpos   = [Hpos_tot(1:3, :); Hpos_tot(5:6, :);Hpos_tot(9, :)]./1E-12;

     Hpos2_tot = computeSecondOrderPartials_FD(n_max, normalized,...
         Cnm_t, Snm_t, Re, GM, rn_ACI, ACAF_ACI, ACAF_ACI);
     Hpos2     = Hpos2_tot./1E-12;
     Hpos2(abs(Hpos2) < 1E-5) = 0;

    % compute Higher order residuals [mE]
    dr_mat = deltaR(:, j) * deltaR(:, j)';
    dr     = [dr_mat(1,1);dr_mat(1,2);dr_mat(1,3);...
              dr_mat(2,2);dr_mat(2,3);dr_mat(3,3)];
    dY_HOT(:, j) = (Y_true(:, j) - y_nom) - (Hpos * deltaR(:, j));
    dY_FOT(:, j) = Hpos * deltaR(:, j);
    dY_SOT(:, j) = Hpos2 * dr;

    % analytical mean
    mean_analy_HOT(:, j) = Hpos2 * s_mean;
end

% compute mean for the HOT
mean_num_HOT = mean(dY_HOT');


figure() 
tt = ["xx", "xy", "xz", "yy", "yz", "zz"];
for k = 1:6
    subplot(2, 3, k);
% %     plot(time, dY_FOT(k, :), 'Marker','*', 'LineStyle', 'none', ...
% %         'Color', ' r');
    hold all;
    plot(time, dY_SOT(k, :), 'Marker','*', 'LineStyle', 'none', ...
    'Color', ' g');
    plot(time, dY_HOT(k, :), 'Marker','*', 'LineStyle', 'none', ...
        'Color', ' b');

    title(tt(k)); grid on;
end
sgtitle('Higher order effects');

disp('STD FOT = ' + string(std(dY_FOT')));
disp('STD SOT = ' + string(std(dY_SOT')));

disp('MEAN HOT = ' + string(mean(dY_HOT')));



%% FUNCTIONS
function [B_comp, B_norm] = ...
    highOrderResidualBoundsUKF(t, planetParams, r_nom, vn, Ps, y_nom, Hk,...
    ACAF_ACI, Cnm_t, Snm_t)
    % eigendecomposition of Ps
    [U, L] = eig(Ps);
    lam = diag(L);
    
    Nk = [];
    c_sigma = 3;

    % build 6 deterministic error samples at the ellipsoid boundary
    n = 3;
    m = length(y_nom);
    eps_mat = zeros(m, 2*n);  % store residuals
    
    idx = 1;
    for i = 1:n
        e_vec = sqrt(lam(i)) * U(:,i);  % 1-sigma principal axis
        for sign = [+1, -1]
            dr = sign * c_sigma * e_vec;         % e.g., 3-sigma along axis i
            r_i = r_nom + dr;
            

            [Y_ACI, ~] = gradiometer_meas(t ,planetParams, ACAF_ACI,...
                     [r_i', vn'], zeros(9,1), Cnm_t, Snm_t);
             y_i = [Y_ACI(1:3); Y_ACI(5:6);Y_ACI(9)]./1E-12; % [mE]
            
            if ~isempty(Nk)
                % robust: remove first-order via NSM
                eps_i = Nk * (y_i - y_nom);
            else
                % explicit first-order subtraction (be careful w/precision)
                eps_i = (y_i - y_nom) - Hk * dr;
            end
            
            eps_mat(:, idx) = eps_i;
            idx = idx + 1;
        end
    end
    
    % component-wise bound: max magnitude over all 6 points
    B_comp = max(abs(eps_mat), [], 2);
    
    % norm-wise bound
    B_norm = max(vecnorm(eps_mat, 2, 1));  % max over columns
end
