clear;
clc;
close all;

addpath("functions/"); addpath("data/");

set(0,'defaultAxesFontSize',16);
%%                  GRAVITY ESTIMATION WITH COLORED NOISE
% Description: Test the gravity estimation process when colored noise is
% introduced.
% Author: Sergio Coll
% Date: 12/04/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% trajectory errors [m]
sigmaP = 90; Nt = length(time);
deltaR = normrnd(0, sigmaP, [3, Nt]);  % trajectory error

% compute S/C orientation
disp('Computing S/C orientation ...')
[BN_mat] = compute_orientation(time, state, instrumentParams);
disp('  DONE ...')

% generate Measurements
disp('Simulating measurements ...')
state_t = state; state_t(:, 1:3) = state(:, 1:3) + deltaR';
[Y_true, ~] = compute_measurements(instrumentParams, planetParams, ...
    poleParams, time, state_t, Cnm_t, Snm_t, BN_mat);

plot_measurements(time, Y_true)
disp('  DONE ...')


%%                      GRAVITY ESTIMATION

% extract asteroid parameteres
GM = planetParams(1); Re = planetParams(2);
normalized = planetParams(4); 
t = time;
n_max = planetParams(3);

% extract asteroid pole parameters
W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
DEC = poleParams(4);

% nominal trajectory
rn = state(:, 1:3)'; vn = state(:, 4:6)';

% measurement noise
sigma_W = instrumentParams(1, 2) * sqrt(instrumentParams(1, 5));

% gravity field apriori
P0 = eye(Ncs-1).*1E-2;
Pc = eye(3).*(sigmaP^2);

s_mean = [Pc(1,1);Pc(1,2);Pc(1,3);Pc(2,2);Pc(2,3);Pc(3,3)];
Css = 1*cov_second_order(Pc);
Pcc = zeros(6*Nt, 6*Nt);
for k = 1:Nt
    r = (k-1)*6 + (1:6);               % row block index
    c = (k-1)*6 + (1:6);               % column block index
    Pcc(r,c) = Css;
end
c  = zeros(6, 1)*sigmaP^2;
s0 = repmat(s_mean, Nt, 1);

% compute weight measurements matrix
sigma_RW = sqrt(Qb(1,1));
[R] = compute_whitening(Nt, sigma_W, sigma_RW);

disp('Starting gravity determination ...')
Mc = 1;
X_final = nan(Ncs, Mc);
for mc = 1:Mc
    disp("  Mc = " + string(mc) + "/" + string(Mc));

    % simulate measurements
    [Y_true, ~] = compute_measurements(instrumentParams, planetParams, ...
        poleParams, time, state_t, Cnm_t, Snm_t, BN_mat);
    
    % compute nominal grav. coefficients
    X0  = [mat2list(Cnm_t, Snm_t, Nc, Ns)];
    for j = 1:Ncs-1
        X0(j+1) = X0(j+1) + normrnd(0, sqrt(P0(j,j)), [1,1]);
    end
    [Cnm, Snm] = list2mat(n_max, Nc, Ns, X0);

    count   = 0;
    iterMax = 5;
    err     = 1;
    xnot_L  = zeros(Ncs-1, 1);
    err_mat = nan(iterMax, 1); pref_HOT = zeros(6, 6, Nt);
    while (count < iterMax) && (err > 1E-5)
        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;

        % measurement vectors
        pref_vec  = nan(6, Nt);
        H_vec     = nan(6, Ncs-1, Nt);
        Hn_vec    = nan(6, 3, Nt); Hn2_vec    = nan(6, 6, Nt);
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt        = W0 + W * t(j);
            ACAF_ACI  =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)]./1E-12;
        
             % compute position and attitude partials
            [Hpos_tot] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos   = [Hpos_tot(1:3, :); Hpos_tot(5:6, :);Hpos_tot(9, :)]./1E-12;

            % second order position partials
            H2 = computeSecondOrderPartials_FD(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos2 = H2./1E-12;

            % prefit residuals
            y_nom = [Y_ACI(1:3); Y_ACI(5:6);Y_ACI(9)]./1E-12;
            dy = Y_true(:, j) - y_nom;
    
            pref_vec(:, j)    = dy;
            H_vec(:, :, j)    = Hc;
            Hn_vec(:, :, j)   = Hpos;
            Hn2_vec(:, :, j)  = Hpos2;
        end

        % apply withening
        [pref_W, H_W, H_HOT] = apply_whitening(R, pref_vec, H_vec, ...
            Hn_vec, Hn2_vec);
        
        Rn   = diag(H_HOT * s0).^2;
        dY   = pref_W;
        Ax_L = Ax_L + (H_W' * (Rn\H_W));
        Nx_L = Nx_L + (H_W' * (Rn\dY));
        Mxc  = H_W' * (Rn\H_HOT);

        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X0(2:end) = X0(2:end) + XNOT_L(1:Ncs-1);
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X0(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L); disp(err)
        err_mat(count + 1) = err;
    
        % update counter
        count = count + 1;
    end
    
    % final state & formal error
    X_final(:, mc) = X0;
    Px = inv(Ax_L);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pcc*Sxc';
end
disp('  DONE');

error        = abs(X_final - X_true);
formal_error = sqrt(diag(Px)); 

% plot results
figure()
semilogy(1:Ncs-1, abs(X_true(2:end)), 'LineWidth', 2, 'Color', 'k');
hold all;
semilogy(1:Ncs-1, error(2:end, :), 'LineStyle', 'none', 'Color', ...
    'r', 'Marker','*');
semilogy(1:Ncs-1, 3.*formal_error, 'LineStyle', '-', ...
    'Color', 'r', 'LineWidth',2);



%% FUNCTIONS
function Css = cov_second_order(Cpp)
% Computes Cov( [x^2, xy, xz, y^2, yz, z^2]' )
% where [x y z]' ~ N(0, Cpp)
%
% INPUT:
%   Cpp : 3x3 covariance of [x y z]'
%
% OUTPUT:
%   Css : 6x6 covariance of the symmetric quadratic terms

Cxx = Cpp(1,1); Cxy = Cpp(1,2); Cxz = Cpp(1,3);
Cyy = Cpp(2,2); Cyz = Cpp(2,3);
Czz = Cpp(3,3);

Css = zeros(6,6);

% ----- variances -----
Css(1,1) = 2*Cxx^2;
Css(2,2) = Cxx*Cyy + Cxy^2;
Css(3,3) = Cxx*Czz + Cxz^2;
Css(4,4) = 2*Cyy^2;
Css(5,5) = Cyy*Czz + Cyz^2;
Css(6,6) = 2*Czz^2;

% ----- cross-covariances -----
Css(1,2) = 2*Cxx*Cxy;
Css(1,3) = 2*Cxx*Cxz;
Css(1,4) = 2*Cxy^2;
Css(1,5) = 2*Cxy*Cxz;
Css(1,6) = 2*Cxz^2;

Css(2,3) = Cxx*Cyz + Cxy*Cxz;
Css(2,4) = 2*Cxy*Cyy;
Css(2,5) = Cxy*Cyz + Cyy*Cxz;
Css(2,6) = Cxz*Cyz;

Css(3,4) = 2*Cxz*Cyy;
Css(3,5) = Cxx*Cyz + Cxy*Czz;
Css(3,6) = 2*Cxz*Czz;

Css(4,5) = 2*Cyy*Cyz;
Css(4,6) = 2*Cyz^2;

Css(5,6) = 2*Cyz*Czz;

% enforce symmetry
Css = Css + triu(Css,1).';

end

