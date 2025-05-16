clear;
clc;
close all;
format long g;

addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_navigation/data/')
addpath('../data/')
addpath('../../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  GRAVITY FIELD DETERMINATION WITH NSM                 %
% Author: Sergio Coll Ibars                                             %
% Date: 05/09/2025                                                      %
% Description: Do gravity field estimation accounting for positon and   %
% attitude errors and use the NSM, LS or a comparison of both           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters
Planet   = "Bennu";       % options: Earth / Bennu / Eros
Solver   = "NSM";        % options: NSM / LS / Both
Errors   = "attitude";    % options: position / attitude / both
saveData = 0;             % options: 0 / 1

[planetParams, poleParams, Kaula, r, Xtrue] = loadPlanet(Planet);
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
DEC = poleParams(4);

% Time options
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 3;
f   = 1/10;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);

% Estimation parameters
P0  = diag(Kaula(2:end).^2);  
S   = normrnd(0, 1E-40 * Kaula);
Xp  = Xtrue + S'; 
[Cnm, Snm] = list2mat(n_max, Nc, Ns, Xtrue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing trajectory ... ')
% Integrate trajectory
if(saveData == 0)
    r0 = [r;0;0];           % [ACI]
    v0 = [0;0;sqrt(GM/r)];  % [ACI]
    
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    PHI0 = reshape(eye(6,6), [36, 1]);
    [~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 4, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
else
    load('Nov_L2position.mat');   % ECEF coordinates
    load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF

    rn_ECEF = positions(:, 2:end)';
    t = positions(:, 1);
    Nt = length(t);

    rn_ACI = rotate2ECI(rn_ECEF, t);
    state_t = zeros(Nt, 6);
    state_t(:, 1:3) = rn_ACI';
end

% compute attitude states
[ACAF_ACI] = compute_planetAtt(poleParams, Planet, saveData , t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position and attitude error
Ar = zeros(3, Nt);

rn = state_t(:, 1:3)' + Ar;
vn = state_t(:, 4:6)';

Amp      = 4.85E-5;
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
disp('Generating measurements ... ')
% mesurement noise & weight
sigma    = 1E-30;                   % [1/s^2]
means    = zeros(1, 9);
std_devs = sigma * ones(1, 9); 
num_realizations = length(t);       % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

R = diag(std_devs.^2);

% compute measurements
[Ytrue, ~, ~] = gradiometer_meas(t ,planetParams, ACAF_ACI, state_t, ...
                zeros(9, Nt), Cnm, Snm, eye(3,3));
[Ytrue] = add_angularComponents(Ytrue, theta, At, flipud(angVel_true),...
    flipud(angAcc_true));  
Ytrue  = Ytrue + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running gravity determination ... ')
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
        [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, R, P0, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
    else
        % TBD
    end
elseif(Solver == "LS")      % run standard LS only
    if(Errors == "attitude")
        [SH_N, sigma_N] = LS_solver_att(planetParams, ACAF_ACI, R, P0, Pc, Pxc, Xp, t ,...
            theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
    else
        % TBD
    end
elseif(Solver == "Both")    % run NSM & LS
    disp(' Solving with NSM .....')
    [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, R, P0, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);

    disp(' Solving with LS .....')
    [SH_LS, sigma_LS] = LS_solver_att(planetParams, ACAF_ACI, R, P0, Pc, Pxc, Xp, t ,...
        theta, thetaDot, thetaDdot, angVel_nom, angAcc_nom, Ytrue, rn, vn);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot resutls
if(Solver == "Both")
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
function [rn_ECI] = rotate2ECI(rn_ECEF, t)
    rn_ECI = rn_ECEF.*0;

    for j = 1:length(t)
        % rotation matrix
        [ECI_ECEF] = convert_ECEF2ECI(t(j));

        rn_ECI(:, j) = ECI_ECEF * rn_ECEF(:, j);
    end
end

function [ACAF_ACI] = compute_planetAtt(poleParams, Planet, saveData , t)
    Nt = length(t);
    ACAF_ACI = ones(3*Nt, 3) * NaN;
    if(Planet == "Earth" && saveData == 1)
        for j =  1:Nt
            [ACI_ACAF] = convert_ECEF2ECI(t(j));
            maxPos = 3*j; minPos = maxPos - 2;
            ACAF_ACI(minPos:maxPos, :) = ACI_ACAF';
        end
    else
        W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
        DEC = poleParams(4);
        for j =  1:Nt
            Wt = W0 + W * t(j);
            R =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            maxPos = 3*j; minPos = maxPos - 2;
            ACAF_ACI(minPos:maxPos, :) = R;
        end
    end
end