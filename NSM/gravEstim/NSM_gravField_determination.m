clear;
clc;
close all;
format long g;

setpref('Internet','E_mail','sergiocollibars@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','sergiocollibars@gmail.com');
setpref('Internet','SMTP_Password','ptcq ybug jcgh wajg');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true'); % For TLS
props.setProperty('mail.smtp.port','587'); % Common TLS port
email = 0;  % send plots by email 1 = yes / 0 = no

% Start diary logging
diaryFile = fullfile(tempdir, 'console_log.txt');
if isfile(diaryFile)
    delete(diaryFile);
end
diary(diaryFile);
diary on;

addpath(genpath('../functions/'))
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
Solver   = "NSM";         % options: NSM / LS / Both
Errors   = "attitude";    % options: position / attitude / both
measMode = "2";           % options 1 = GG + ST + GYRO / 2 = GG + ST 
saveData = 0;             % options: 0 / 1

type_Attitude = "RTN";    % options: inertial / RTN / GRF
type_posErr = "constant"; % options: constant / periodic / random
type_attErr = "random";   % options: constant / periodic / linear

printConfig_console(Planet, Solver, Errors, measMode, saveData, ...
    type_Attitude, type_posErr, type_attErr);

%                Mesurement mask 
%           xx xy xz yx yy yz zx zy zz
mask     =  [1, 1, 1, 1, 1, 1, 1, 1, 1]';
printMask_console(mask);

[planetParams, poleParams, Kaula, r, Xtrue, Iner] = loadPlanet(Planet);
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
DEC = poleParams(4);

% Time options
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 20 * 30;
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
[orientation, Mext] = SC_orientation(t, state_t, type_Attitude,...
    GRF_ACI, Iner);

% attitude nominal value
theta    = orientation(1:3, :);  % Euler angle      [rad]
thetaDot = orientation(4:6, :);  % Euler angle rate [rad /s]

% compute Body to ACI rotation matrix
B_ACI_mat = ones(3*Nt, 3);
for j = 1:length(t)
    maxVal = 3 * j; minVal = maxVal -2;
    
    yaw  = theta(1, j); pitch = theta(2,j); roll = theta(3, j);
    B_ACI =rotationMatrix(yaw, pitch, roll, ...
        [3, 2, 1]);

    B_ACI_mat(minVal:maxVal, :) = B_ACI;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating state errors ... ')
% Position error
Amp  = 0.*[1;0.7;0.5];          % [m] 
[Ar] = generate_posErrors(t, type_posErr, Amp, T);

% Attitude error
Amp  = 1.45E-5.*[1;1;1];        % [rad] 
[att_Err] = generate_attErrors(t, type_attErr, Planet, saveData, Amp, T, measMode);

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
% measurement noise & weight
sigma    = 1E-12;                   % [1/s^2]
means    = zeros(1, 9);
std_devs = sigma * ones(1, 9); 
num_realizations = length(t);       % Number of realizations

noise = normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));
% % load("noise.mat");

R = diag(std_devs(logical(mask)).^2);

% compute measurements
[Ytrue, ~, ~] = gradiometer_meas(t ,planetParams, ACAF_ACI, state_t, ...
                zeros(9, Nt), Cnm, Snm);
[Ytrue] = add_angularComponents(Ytrue, theta, att_Err(1:3, :), angVel_true,...
    angAcc_true); 
Ytrue  = Ytrue(logical(mask), :) + noise(logical(mask), :);

% plot gradiometer signal
if(length(Ytrue(:, 1)) > 6)
    plot_signal(Ytrue, t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running gravity determination ... ')
% filter for NaN values in position ,velocity or attitude vector
A = [rn;vn; theta; angVel_nom; angAcc_nom];
nanMask = isnan(A);               % Logical matrix: true where A is NaN
colHasNaN = any(nanMask, 1);      % 1 means operate across rows → result is 1 x N

% Now find the first column where no NaNs exist
firstIndx = find(~colHasNaN, 1, 'first');

if(Errors == "attitude")
    err = [att_Err(1:3, firstIndx:end);...
        angVel_true(:, firstIndx:end) - angVel_nom(:, firstIndx:end)];
    Pc = zeros(6, 6); Pxc = zeros(Ncs - 1, 6);
    for j = 1:3
        Pc(j, j) = err(j, 1)^2;
    end
    for j = 4:6
        Pc(j, j) = err(j, 1)^2;
    end
else
    % TBD
end

% short the input values to remove NaNs
maxPos = 3 * firstIndx; minPos = maxPos -2;
ACAF_ACI    = ACAF_ACI(minPos:end, :);
B_ACI_mat   = B_ACI_mat(minPos:end, :);
t           = t(1, firstIndx:end);
theta       = theta(:, firstIndx:end);
angVel_nom  = angVel_nom(:, firstIndx:end);
angAcc_nom  = angAcc_nom(:, firstIndx:end);
Ytrue       = Ytrue(:, firstIndx:end);
rn          = rn(:, firstIndx:end);
vn          = vn(:, firstIndx:end);

% run estimation
if(Solver == "NSM")         % run NSM only
    if(Errors == "attitude")
        [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, B_ACI_mat,...
            R, P0, Pc, Pxc, Xp, t , theta, angVel_nom, angAcc_nom, Ytrue, rn,...
            vn, Iner, mask);
    else
        % TBD
    end
elseif(Solver == "LS")      % run standard LS only
    if(Errors == "attitude")
        [SH_N, sigma_N] = LS_solver_att(planetParams, ACAF_ACI, B_ACI_mat,...
            R, P0, Pc, Pxc, Xp, t ,...
            theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner, mask);
    else
        % TBD
    end
elseif(Solver == "Both")    % run NSM & LS
    disp('Solving with NSM .....')
    [SH_N, sigma_N] = NSM_solver_att(planetParams, ACAF_ACI, B_ACI_mat,...
        R, P0, Pc, Pxc, Xp, t ,...
        theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner, mask);

    disp('Solving with LS .....')
    [SH_LS, sigma_LS] = LS_solver_att(planetParams, ACAF_ACI, B_ACI_mat,...
        R, P0, Pc, Pxc, Xp, t ,...
        theta, angVel_nom, angAcc_nom, Ytrue, rn, vn, Iner, mask);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Displaying results ...')
% plot resutls
if(Solver == "Both")
    if(n_max > 10)
        limit = [1E-2, 50];
        limitTicks = [1E-2, 1E-1, 1, 5, 20, 50];
        scale = 'log'; digits = '%.3f';

        % plot gravity field error pr SH coefficient. Pyramide plot
        tt = "NSM estimation error "; 
        error = (abs(Xtrue - [1;SH_N])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt, limit, ...
            limitTicks, scale, digits);

        tt = "LS estimation error "; 
        error = (abs(Xtrue - [1;SH_LS])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt, limit,...
            limitTicks, scale, digits);
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
        limit = [1E-2, 50];
        limitTicks = [1E-2, 1E-1, 1, 5, 20, 50];
        scale = 'log'; digits = '%.3f';

        % plot gravity field error pr SH coefficient. Pyramide plot
        tt = Solver +  " estimation error "; 
        error = (abs(Xtrue - [1;SH_N])./abs(Xtrue)).* 100;
        plot_high_gravField(n_max, Nc, Ns, error, tt, limit,...
            limitTicks, scale, digits);
    else
        mk = 'square';  tt = "Estimation error "; 
        llg = {'truth','NSM', 'error'};
        plot_gravField(Xtrue, sigma_N, Xtrue(2:end) - SH_N, n_max, tt, mk, llg);
    end

    % plot gravity field error per RMS value
    tt = "RMS error using NSM"; llg = {'truth', '3 \sigma NSM', 'NSM error'};
    error = Xtrue - [1;SH_N];
    plot_gravField_RMS(Xtrue, sigma_N, error, n_max, tt, llg);
end

diary off; % close diary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(email)
        % load console output 
        consoleText = fileread(diaryFile);

        sendOpenFiguresByEmail(...
            'sergiocollibars@gmail.com', ...
            'Auto Report - MATLAB Figures', consoleText);
end
delete(diaryFile);  % clean up immediately