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
Planet   = "Earth";       % options: Earth / Bennu / Eros
Solver   = "Both";        % options: NSM / LS / Both
Errors   = "attitude";    % options: position / attitude / both
measMode = "2";           % options 1 = GG + ST + GYRO / 2 = GG + ST 
saveData = 0;             % options: 0 / 1

[planetParams, poleParams, Kaula, r, Xtrue, Iner] = loadPlanet(Planet);
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
DEC = poleParams(4);

% Time options
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 10;
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
    ACI_ECEF = 0; GRF_ACI = 0;
else
    load('Nov_L2position.mat');   % ECEF coordinates
    load('Nov_L2velocity.mat');   % ECEF coordinates
    load('Nov_L2ECEF2ITRF.mat');  % rotation matrix ECEF 2 ITRF
    load('Nov_L2ITRF2GRF.mat');   % rotation matrix ITRF 2 GRF
    
    [~, t2]  = quaternion2CDM(outPut);
    [~, t3]  = quaternion2CDM(outPut2);
    
    % Check time-points to make sure that all files are at the
    % same time.
    t1 = positions(:, 1);
    commonTimes = intersect(intersect(t1, t2), t3);
    
    [~, idx1] = ismember(commonTimes, t1);
    [~, idx2] = ismember(commonTimes, t2);
    [~, idx3] = ismember(commonTimes, t3);

    t = t1(idx1)';
    Nt = length(t);
    
% %     % WARNING: reducing data set
% %     t = t(1:round(Nt/10));
% %     Nt = length(t);

    rn_ECEF = positions(idx1, 2:end)';
    vn_ECEF = velocity(idx1, 2:end)';
    ACI_ECEF = quaternion2CDM(outPut(idx2, :));
    GRF_ACI  = quaternion2CDM(outPut2(idx3, :));

    state_ACI = rotate2ECI([rn_ECEF; vn_ECEF], ACI_ECEF, t);
    state_t(:, 1:3) = state_ACI(1:3, :)';
    state_t(:, 4:6) = state_ACI(4:6, :)';
end

% compute planet orientation parameters
[ACAF_ACI] = compute_planetAtt(poleParams, t, saveData, ACI_ECEF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Computing attitude ... ')
type = "GRF";
[orientation, Mext] = SC_orientation(t, state_t, type, GRF_ACI, Iner);

% check for outliers
% % [orientation(7:9, :)] = check_outliers(t, orientation(7:9, :));

% attitude nominal value
theta    = orientation(1:3, :);  % Euler angle      [rad]
thetaDot = orientation(4:6, :);  % Euler angle rate [rad /s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating state errors ... ')
% Position error
type = "constant";              % options: constant / periodic / random
Amp  = 0.*[1;0.7;0.5];          % [m] 
[Ar] = generate_posErrors(t, type, Amp, T);

% Attitude error
type = "FOGMP";              % options: constant / periodic / linear
Amp  = 4.85E-10.*[1;0.7;0.5];  % [rad] 
[att_Err] = generate_attErrors(t, type, Planet, saveData, Amp, T, measMode);

% include position errors
rn = state_t(:, 1:3)' + Ar;
vn = state_t(:, 4:6)';

% include attitude errors
[angVel_true, angAcc_true] = compute_angularVals(theta, thetaDot, Iner,...
    measMode, att_Err, Mext);
[angVel_nom, angAcc_nom]   = compute_angularVals(theta, thetaDot, Iner,...
    measMode, zeros(6, Nt), Mext);

% plot S/C state
err_omega    = angVel_true - angVel_nom;
err_omegaDot = angAcc_true - angAcc_nom;
plot_SC_state(t, saveData, Planet, state_t, orientation, ...
    angVel_nom, angAcc_nom, err_omega, err_omegaDot);

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
    err = [att_Err(1:3, :);angVel_true - angVel_nom];
    Pc = zeros(6, 6); Pxc = zeros(Ncs - 1, 6);
    for j = 1:3
        Pc(j, j) = max(err(j, :))^2;
    end
    for j = 4:6
        Pc(j, j) = max(err(j-3, :)*1E-3)^2;
    end
else
    % TBD
end

% run estimation
if(Solver == "NSM")         % run NSM only
    if(Errors == "attitude")
        [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, R, P0, Xp, t ,...
        theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner);
    else
        % TBD
    end
elseif(Solver == "LS")      % run standard LS only
    if(Errors == "attitude")
        [SH_N, sigma_N] = LS_solver_att(planetParams, ACAF_ACI, R, P0, Pc, Pxc, Xp, t ,...
            theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner);
    else
        % TBD
    end
elseif(Solver == "Both")    % run NSM & LS
    disp(' Solving with NSM .....')
    [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, R, P0, Xp, t ,...
        theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner);

    disp(' Solving with LS .....')
    [SH_LS, sigma_LS] = LS_solver_att(planetParams, ACAF_ACI, R, P0, Pc, Pxc, Xp, t ,...
        theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Displaying results ...')
% plot resutls
if(Solver == "Both")
    if(n_max > 10)
        % plot gravity field error pr SH coefficient. Pyramide plot
        tt = "NSM estimation error "; 
        error = (abs(Xtrue - [1;SH_N])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt);

        tt = "LS estimation error "; 
        error = (abs(Xtrue - [1;SH_LS])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt);
    else
        % plot gravity field error per SH coefficient
        mk = 'square';  tt = "NSM estimation error "; 
        llg = {'truth','NSM', 'error'};
        plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);
    
        mk = 'square';  tt = "LS estimation error "; 
        llg = {'truth','LS', 'error'};
        plot_gravField(Xtrue, sigma_LS, Xtrue(2:end) - SH_LS, n_max, tt, mk, llg);
    end

    % plot gravity field error per RMS value
    tt = "RMS error using NSM"; llg = {'truth', '3 \sigma NSM', 'NSM error'};
    error = Xtrue - [1;SH_N];
    plot_gravField_RMS(Xtrue, sigma_N, error, n_max, tt, llg);

    tt = "RMS error using LS"; llg = {'truth', '3 \sigma LS', 'LS error'};
    error = Xtrue - [1;SH_LS];
    plot_gravField_RMS(Xtrue, sigma_LS, error, n_max, tt, llg);

else
    if(n_max > 10)
        % plot gravity field error pr SH coefficient. Pyramide plot
        tt = Solver +  " estimation error "; 
        error = (abs(Xtrue - [1;SH_N])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt);
    else
        mk = 'square';  tt = "Estimation error "; 
        llg = {'truth','NSM', 'error'};
        plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
function [output] = check_outliers(t, input)
    ths = 5E-7; 
    for j = 1:length(t)
        if(rms(input(:, j)) > ths)
            input(:, j) = zeros(3, 1);
        end
    end
    output = input;
end