clear;
clc;
close all;
format long g;

addpath('../functions/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_navigation/data/')
addpath('../data/')
addpath('../../matlab_codes/GOCE_products/GOCE_L2b_MatlabReaders/data/')
addpath('../../QGG_gravEstim/data_files/')
addpath('../../QGG_gravEstim/src/')

set(0,'defaultAxesFontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  GRAVITY FIELD DETERMINATION WITH NSM                 %
% Author: Sergio Coll Ibars                                             %
% Date: 05/09/2025                                                      %
% Description: Do gravity field estimation accounting for positon and   %
% attitude errors and use the NSM, LS or a comparison of both           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters
Planet   = "Bennu";       % options: Earth / Bennu / Eros
Solver   = "Both";        % options: NSM / LS / Both
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
S   = normrnd(0, 1E-1 * Kaula);
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
    ACI_ECEF = 0;
else
    load('Nov_L2position.mat');   % ECEF coordinates
    load('Nov_L2velocity.mat');   % ECEF coordinates
    load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF
    
    [ACI_ECEF] = read_ECEF2ITRF_mat(outPut);

    rn_ECEF = positions(:, 2:end)';
    vn_ECEF = velocity(:, 2:end)';
    t = positions(:, 1)';
    Nt = length(t);

    rn_ACI = rotate2ECI(rn_ECEF, ACI_ECEF, t);
    vn_ACI = rotate2ECI(vn_ECEF, ACI_ECEF, t);
    state_t = zeros(Nt, 6);
    state_t(:, 1:3) = rn_ACI';
    state_t(:, 4:6) = vn_ACI';
end

% compute planet orientation parameters
[ACAF_ACI] = compute_planetAtt(poleParams, t, saveData, ACI_ECEF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing attitude ... ')
type = "RTN";
[orientation] = SC_orientation(t, state_t, type);

% check for outliers
% % [orientation(7:9, :)] = check_outliers(t, orientation(7:9, :));

% attitude nominal value
theta    = orientation(1:3, :);  % [rad]
thetaDot = orientation(4:6, :);  % [rad /s]
thetaDdot= orientation(7:9, :);  % [rad / s^2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating state errors ... ')
% Position error
type = "constant";              % options: constant / periodic / random
Amp  = 0.*[1;0.7;0.5];          % [m] 
[Ar] = generate_posErrors(t, type, Amp, T);

% Attitude error
type = "periodic";              % options: constant / periodic / linear
Amp  = 4.85E-10.*[1;0.7;0.5];    % [rad] 
[att_Err] = generate_attErrors(t, type, Planet, Amp, T);

% include position errors
rn = state_t(:, 1:3)' + Ar;
vn = state_t(:, 4:6)';

% include attitude errors
[angVel_true, angAcc_true] = compute_angularVals(theta + att_Err(1:3, :), ...
    thetaDot + att_Err(4:6, :), thetaDdot + att_Err(7:9, :));
[angVel_nom, angAcc_nom]   = compute_angularVals(theta, thetaDot, thetaDdot);

% plot S/C state
plot_SC_state(t, Planet, state_t, orientation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating measurements ... ')
% mesurement noise & weight
sigma    = 1E-12;                   % [1/s^2]
means    = zeros(1, 9);
std_devs = sigma * ones(1, 9); 
num_realizations = length(t);       % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

R = diag(std_devs.^2);

% compute measurements
[Ytrue, ~, ~] = gradiometer_meas(t ,planetParams, ACAF_ACI, state_t, ...
                zeros(9, Nt), Cnm, Snm, eye(3,3));
[Ytrue] = add_angularComponents(Ytrue, theta, att_Err(1:3, :), angVel_true,...
    angAcc_true);  
Ytrue  = Ytrue + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running gravity determination ... ')
if(Errors == "attitude")
    Pc = zeros(6, 6); Pxc = zeros(Ncs - 1, 6);
    for j = 1:3
        Pc(j, j) = max(att_Err(j, :))^2;
    end
    for j = 4:6
        Pc(j, j) = max(att_Err(j-3, :)*1E-3)^2;
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
disp('Displaying results ...')
% plot resutls
if(Solver == "Both")
    % plot gravity field error per SH coefficient
    mk = 'square';  tt = "NSM estimation error "; llg = {'truth','NSM', 'error'};
    plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);

    mk = 'square';  tt = "LS estimation error "; llg = {'truth','LS', 'error'};
    plot_gravField(Xtrue, sigma_LS, Xtrue(2:end) - SH_LS, n_max, tt, mk, llg);

    % plot gravity field error per RMS value
    tt = "RMS error using NSM"; llg = {'truth', '3 \sigma NSM', 'NSM error'};
    plot_gravField_RMS(Xtrue, sigma_N, SH_N, n_max, tt, llg);

    tt = "RMS error using LS"; llg = {'truth', '3 \sigma LS', 'LS error'};
    plot_gravField_RMS(Xtrue, sigma_LS, SH_LS, n_max, tt, llg);

else
    mk = 'square';  tt = "Estimation error "; llg = {'truth','NSM', 'error'};
    plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
function [ACAF_ACI] = compute_planetAtt(poleParams , t, saveData, ACI_ECEF)
    Nt = length(t);
    ACAF_ACI = ones(3*Nt, 3) * NaN;
    if(saveData == 1)
        for j = 1:Nt
            maxPos = 3 * j; minPos = maxPos - 2;
            ACAF_ACI(minPos:maxPos, :) = ACI_ECEF(minPos:maxPos, :)';
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

function [output] = check_outliers(t, input)
    ths = 5E-7; 
    for j = 1:length(t)
        if(rms(input(:, j)) > ths)
            input(:, j) = zeros(3, 1);
        end
    end
    output = input;
end