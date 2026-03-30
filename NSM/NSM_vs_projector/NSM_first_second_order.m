clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Test the NSM and ANSM ability to eliminte the first & second
% order effects.
% Author: Sergio Coll
% Date: 11/16/24

% traget body
target_body = "Bennu";  % options: Bennu / Eros

if(target_body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
    W = 4.06130329511851E-4;          % Rotation ang. vel   [rad/s]
    W0 = 0;                           % Initial asteroid longitude
    RA = deg2rad(86.6388);            % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);          % Declination         [rad]
    r = 0.36E3;                       % orbit radius        [m]
    sigmaP = 0.1;                     % orbit uncertainty   [m]
    sigmaE = 10 * pi / (3600 * 180);  % attitude errors     [rads]
    sigmaW = 1E-3;                    % attitude errors     [deg/sqrt(hr)]
elseif(target_body == "Eros")
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    n_max  = 0;
    normalized = 1;
    GM =  459604.431484721;          % Point mass value    [m^3/s^2]
    W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
    W0 = 0;                          % Initial asteroid longitude
    RA = deg2rad(11.363);            % Right Ascension     [rad]
    DEC = deg2rad(17.232);           % Declination         [rad]
    r = 24E3;                        % orbit radius        [m]
    sigmaP = 3;                      % orbit uncertainty   [m]
    sigmaE = 10 * pi / (3600 * 180); % attitude errors     [rads]
    sigmaW = 1E-3;                   % attitude errors     [deg/sqrt(hr)]
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
f = 1/10;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);
At = t(2) - t(1);                           % [sec]

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

sigmaP_array = linspace(sigmaP,100*sigmaP,3); dX = ones(Ncs-1, 1).*1000;
RMS_res      = ones(6, length(sigmaP_array)) * NaN;
RMS_res_NSM  = ones(3, length(sigmaP_array)) * NaN; 
RMS_res_ANSM = ones(2, length(sigmaP_array)) * NaN;
for k = 1:length(sigmaP_array)
    disp('Error iterataion ' + string(k))

    % state errors
    deltaPos      =  [1;1;1].*sigmaP_array(k).*ones(3, Nt);   % [m]
    
    % output variables
    Res     = ones(6, Nt) * NaN; 
    Res_NSM = ones(3, Nt) * NaN; Res_ANSM = ones(2, Nt) * NaN;
    for j = 1:Nt
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % Inertially-fixed
        rn_ACI = rn(:, j);
        ACAF_BODY = ACAF_ACI;
        B_ACI     = eye(3,3);
    
        w = [0;0;0];
        
        % compute gradiometer measurements
        [Y_ACI, H_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
            [rn(:, j)', vn(:, j)'], zeros(9,1), Cnm, Snm);
    
        % Rotate measurements
        H = [H_ACI(1:3, 2:end); H_ACI(5:6, 2:end); H_ACI(9, 2:end)];
    
        % compute 1st order position partials
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
            rn_ACI, ACAF_ACI, ACAF_BODY);
        Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];
    
        % compute 2nd order position partials
        Hpos_2 = compute_posPartials_2ndOrder(GM, rn_ACI(1), ...
            rn_ACI(2), rn_ACI(3));
        Hpos_2_z     = -Hpos_2(1, :) - Hpos_2(4, :);
        Hpos_2_total = [Hpos_2;Hpos_2_z];
    
        Hpos_total = [Hpos, Hpos_2_total];
        
        % residuals (1st and 2nd order)
        Res_1 = Hpos * deltaPos(:, j);
        dP = deltaPos(:, j) * deltaPos(:, j)';
        Res_2 = Hpos_2_total * [dP(1,1);dP(1,2);dP(1,3);dP(2,2);...
            dP(2,3); dP(3,3)];
    
        Res(:, j) = Res_1 + Res_2;
    
        % compute NSM & residuals
        C_NSM = null(Hpos');
        Res_NSM(:, j) = C_NSM' * Res(:, j);
        H_NSM(:, j) = C_NSM' * H * dX;
    
        % compute ANSM & residuals
        [s,v,d] = svd((C_NSM' * Hpos_2_total)');
        C_ANSM = d(:, 3);
        Res_ANSM(:, j) = C_ANSM' * Res_NSM(:, j);
        H_ANSM(:, j) = C_ANSM' * H_NSM(:, j);
    end

    % compute RMS values
    RMS_res(:, k)      = rms(Res, 2);
    RMS_res_NSM(:, k)  = rms(Res_NSM, 2);
    RMS_res_ANSM(:, k) = rms(Res_ANSM, 2);
end

figure()
plot(sigmaP_array, mean(RMS_res), 'LineWidth', 2);
xlabel('Position error \sigma_P')
ylabel('Residuals RMS');
title('Position error residual, original space')

figure()
plot(sigmaP_array, mean(RMS_res_NSM),  'LineWidth', 2);
hold all;
plot(sigmaP_array, RMS_res_ANSM, 'LineWidth', 2);
grid on;
xlabel('Position error \sigma_P')
ylabel('Residuals RMS');
title('Position error residual, NSM & ANSM space')
legend('NSM', 'ANSM')