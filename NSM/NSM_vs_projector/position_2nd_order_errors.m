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
    n_max  = 10;
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

% state errors
sigmaW        = sigmaW / sqrt(At) * pi / (180*3600);                        % [rad/s]
deltaPos      = normrnd(0, sigmaP, [3, Nt]) + [1;1;1].*0.01;                % [m]
deltaAtt      = normrnd(0, sigmaE, [3, Nt]) + [1;6;10].*(pi/(180*3600));    % [rad]
deltaOmega    = normrnd(0, sigmaW, [3, Nt]);                                % [rad/s]
deltaOmegaDot = (deltaOmega./At);                                           % [rad/s^2]

% S/C inertia matrix
Iner = [580, 0, 0;...
    0,640,0;...
    0,0,1100];

% output variables
Res_1_pos = ones(6, Nt) * NaN; Res_1_att = ones(6, Nt) * NaN;
Res_2_pos = ones(6, Nt) * NaN; Res_2_att = ones(6, Nt) * NaN;
for j = 1:Nt
    disp(j/Nt * 100)
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % Inertially-fixed
    rn_ACI = rn(:, j);
    ACAF_BODY = ACAF_ACI;
    B_ACI     = eye(3,3);

    w = [0;0;0];
    
    % compute gradiometer measurements
    [Y_ACI, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
        [rn(:, j)', vn(:, j)'], zeros(9,1), Cnm, Snm);

    % Rotate measurements
    T_ACI = reshape(Y_ACI, [3,3]);

    % Rotation matrix error
    BN = rotationMatrix(deltaAtt(1, j), deltaAtt(2, j), deltaAtt(3, j), ...
        [3, 2, 1]);

    T_B = BN * T_ACI * BN';
    Y_B = reshape(T_B, [9, 1]);

    % compute 1st order attitude errors
    [Hrot_grad] = compute_rotPartials_analy(Y_ACI, B_ACI);
    [Hrot_omega_dyad, H_omegaDot_dyad, ~, ~] = ...
        compute_angularDyadPartials_v2(w, Iner);

     res = Hrot_grad * deltaAtt(:, j);

     Res_1_att(:, j) = [res(1:3); res(5:6); res(9)];

     % compute 2nd order attitude errors
     [omega_dyad, omegaDot_dyad] = compute_Skew_mat(deltaOmega(:, j), ...
         deltaOmegaDot(:, j));
     res2 = (Y_B - Y_ACI) -  res - reshape(omega_dyad, [9, 1]);
     Res_2_att(:, j) = [res2(1:3); res2(5:6); res2(9)];

    % compute 1st order position partials
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
        rn_ACI, ACAF_ACI, ACAF_BODY);
    Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

    % compute 2nd order position partials
    Hpos_2 = compute_posPartials_2ndOrder(GM, rn_ACI(1), ...
        rn_ACI(2), rn_ACI(3));
    Hpos_2_z     = -Hpos_2(1, :) - Hpos_2(4, :);
    Hpos_2_total = [Hpos_2;Hpos_2_z];
    
    % residuals (1st and 2nd order)
    Res_1_pos(:, j) = Hpos * deltaPos(:, j);
    dP = deltaPos(:, j) * deltaPos(:, j)';
    Res_2_pos(:, j) = Hpos_2_total * [dP(1,1);dP(1,2);dP(1,3);dP(2,2);...
        dP(2,3);dP(3,3)];
end

components = {'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}';

figure()
for k = 1:6
    subplot(2, 3, k)
    plot(t./T, Res_1_pos(k, :)./1E-12, '.', 'Color', 'b')
    hold on;
    plot(t./T, Res_1_att(k, :)./1E-12, '.', 'Color', 'r')
    ylabel('milli-Eotvos')
    grid on;
    xlabel('Orbit rev.')
    title(components(k))
end
sgtitle('1st order residuals')
legend('Pos.', 'Att.')

figure()
for k = 1:6
    subplot(2, 3, k)
    plot(t./T, Res_2_pos(k, :)./1E-12, '.', 'Color', 'b')
    hold on;
    plot(t./T, Res_2_att(k, :)./1E-12, '.', 'Color', 'r')
    ylabel('milli-Eotvos')
    grid on;
    xlabel('Orbit rev.')
    title(components(k));
end
sgtitle('2nd order residuals')
legend('Pos.', 'Att.');

minVals  = min(Res_1_pos')'./1E-12; 
maxVals  = max(Res_1_pos')'./1E-12;
rmsVals  = rms(Res_1_pos, 2)./1E-12;
table(components, minVals, maxVals, rmsVals)


minVals  = min(Res_1_att')'./1E-12; 
maxVals  = max(Res_1_att')'./1E-12;
rmsVals  = rms(Res_1_att, 2)./1E-12;
table(components, minVals, maxVals, rmsVals)


minVals  = min(Res_2_pos')'./1E-12; 
maxVals  = max(Res_2_pos')'./1E-12;
rmsVals  = rms(Res_2_pos, 2)./1E-12;
table(components, minVals, maxVals, rmsVals)


minVals  = min(Res_2_att')'./1E-12; 
maxVals  = max(Res_2_att')'./1E-12;
rmsVals  = rms(Res_2_att, 2)./1E-12;
table(components, minVals, maxVals, rmsVals)


%% FUNCTIONS

function [omega_dyad, omegaDot_dyad] = compute_Skew_mat(angVel, angAcc)
    j = 1;
  omega = [0, angVel(3, j), -angVel(2, j);...
    -angVel(3, j), 0 ,angVel(1, j);...
     angVel(2, j), -angVel(1, j), 0];
  omegaDot = [0, angAcc(3, j), -angAcc(2, j);...
    -angAcc(3, j), 0 ,angAcc(1, j);...
    angAcc(2, j), -angAcc(1, j), 0];
  omega_dyad = omega^2;
  omegaDot_dyad = omegaDot;
end