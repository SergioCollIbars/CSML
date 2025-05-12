clear;
clc;
close all;
format long g;

addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_navigation/data/')
addpath('../data/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  GRAVITY FIELD DETERMINATION WITH NSM                 %
% Author: Sergio Coll Ibars                                             %
% Date: 05/09/2025                                                      %
% Description: Do gravity field estimation accounting for positon and   %
% attitude errors and use the NSM, LS or a comparison of both           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters
Planet = "Earth";       % options: Earth / Bennu / Eros
Solver = "Both";         % options: NSM / LS / Both
Errors = "attitude";    % options: position / attitude / both

[planetParams, poleParams, Kaula, r, Xtrue] = loadPlanet(Planet);
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
DEC = poleParams(4);

% Time options
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 10*16;
f   = 1/10;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);

% Estimation parameters
P0  = diag(Kaula(2:end).^2);  
S   = normrnd(0, 1E-2 * Kaula);
Xp  = Xtrue + S'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Integrate trajectory
r0 = [r;0;0];           % [ACI]
v0 = [0;0;sqrt(GM/r)];  % [ACI]
[Cnm, Snm] = list2mat(n_max, Nc, Ns, Xtrue);

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 4, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position and attitude error
Ar = zeros(3, Nt);

rn = state_t(:, 1:3)' + Ar;
vn = state_t(:, 4:6)';

Amp      = 5E-6;
At       = Amp*ones(3, Nt).*[1;1;1];             % [rad]
dAt_dt   = Amp.*ones(3, Nt).*[0;0;0];            % [rad/s]
ddAt_ddt = zeros(3, Nt);                         % [rad/s^2]

% attitude nominal value
theta    = zeros(3, Nt);                         % nominal attitude [rad]
thetaDot = zeros(3, Nt);
thetaDdot= zeros(3, Nt); 

[angVel_true, angAcc_true] = compute_angularVals(theta + At, ...
    thetaDot + dAt_dt, thetaDdot + ddAt_ddt);
[angVel_nom, angAcc_nom]   = compute_angularVals(theta, thetaDot, thetaDdot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mesurement noise & weight
sigma    = 1E-12;                   % [1/s^2]
means    = zeros(1, 9);
std_devs = sigma * ones(1, 9); 
num_realizations = length(t);       % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

R = diag(std_devs.^2);

% compute measurements
[Ytrue, ~, ~] = gradiometer_meas(t ,planetParams, poleParams, state_t, ...
                zeros(9, Nt), Cnm, Snm, eye(3,3));
[Ytrue] = add_angularComponents(Ytrue, theta, At, flipud(angVel_true),...
    flipud(angAcc_true));  
Ytrue  = Ytrue + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Errors == "attitude")
    Pc = zeros(6, 6); Pxc = zeros(Ncs - 1, 6);
    for j = 1:3
        Pc(j, j) = max(At(j, :))^2;
    end
    for j = 4:6
        Pc(j, j) = max(At(j-3, :))^2;
    end
else
    % TBD
end

% run estimation
if(Solver == "NSM")         % run NSM only
    if(Errors == "attitude")
        [SH_N, sigma_N] = NSM_solver_att(planetParams, poleParams, R, P0, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
    else
        % TBD
    end
elseif(Solver == "LS")      % run standard LS only
    if(Errors == "attitude")
        [SH_N, sigma_N] = LS_solver_att(planetParams, poleParams, R, P0, Pc, Pxc, Xp, t ,...
            theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
    else
        % TBD
    end
elseif(Solver == "Both")    % run NSM & LS
    disp(' Solving with NSM .....')
    [SH_N, sigma_N] = NSM_solver_att(planetParams, poleParams, R, P0, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);

    disp(' Solving with LS .....')
    [SH_LS, sigma_LS] = LS_solver_att(planetParams, poleParams, R, P0, Pc, Pxc, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute RMS values
X_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            Xtrue, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
err_RMS_N  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            Xtrue - [1;SH_N], zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
err_RMS_LS = computeRMS_coeffErr(n_max, Nc, Ns, ...
            Xtrue - [1;SH_LS], zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
s_RMS_N  = computeRMS_coeffErr(n_max, Nc, Ns, ...
           [1;3*sigma_N], zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
s_RMS_LS = computeRMS_coeffErr(n_max, Nc, Ns, ...
            [1;3*sigma_LS], zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 

% plot resutls
if(Solver == "Both")
    mk = 'square';  tt = "NSM estimation error "; llg = {'truth','NSM', 'error'};
    plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);

    mk = 'square';  tt = "LS estimation error "; llg = {'truth','NSM', 'error'};
    plot_gravField(Xtrue, sigma_LS, Xtrue(2:end) - SH_LS, n_max, tt, mk, llg);

    figure()
    semilogy(1:n_max, X_RMS, 'LineWidth', 2, 'Color', 'k');
    hold all;
    semilogy(1:n_max, s_RMS_N, 'LineWidth', 2, 'Color', 'b');
    semilogy(1:n_max, s_RMS_LS, 'LineWidth', 2, 'Color', 'g');
    semilogy(1:n_max, err_RMS_N, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--');
    semilogy(1:n_max, err_RMS_LS, 'LineWidth', 2, 'Color', 'g', 'LineStyle','--');
    legend('truth', '3\sigma NSM', '3\sigma LS', 'NSM error', 'LS error')
    title('RMS error')
    xlim([2, n_max]);

else
    mk = 'square';  tt = "Estimation error "; llg = {'truth','NSM', 'error'};
    plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);
end