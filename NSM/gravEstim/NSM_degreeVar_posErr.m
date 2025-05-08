clear;
clc;
close all;

addpath('../functions/')
addpath('../../QGG_gravEstim/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_gravEstim/data_files/')
addpath('../../QGG_navigation/data/')
set(0,'defaultAxesFontSize',16);

%%            NSM POSITION ERROR EFFECT IN GRAV. ESTIMATION
% Description: Compute the trucation error for position errors depending on
% the order.
% Author: Sergio Coll
% Date: 10/09/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 0;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % % Earth parameters
% % savedData = 0;                % use saved data. 1 = yes / 0 = no
% % path = "HARMCOEFS_EARTH_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % path = "SIGMACOEFS_EARTH_1.txt";
% % [sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
% % GM = 3.986004418E14;
% % n_max  = 20;
% % normalized = 1;
% % W = 2 * pi / (24*3600);     % Rotation ang. vel   [rad/s]
% % W0 = 0;                     % Initial asteroid longitude
% % RA = -pi/2;                 % Right Ascension     [rad]
% % DEC = pi/2;                 % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

% compute Kaula rule
K = 0.026;
[CS_K] = compute_Kaula(n_max, K);

% Initial conditions
r      = 0.4E3;
% % r      = 250E3 + Re; 
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
rev = 6;
f = 1/10;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% weight matrix
noise0 = zeros(9, Nt);
sigmam  = 1E-15;
sigmar = 1;
R  = diag([sigmam, sigmam].^2);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 0, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

% % [S] = mat2list(sigma_Cnm, sigma_Snm, Nc, Ns);
% % S = 1.*S;
% % P0 = diag((S(1:end)).^2);

% Gravity estimation
Ax = 0; Mxc = 0;
for j = 1:Nt   
    disp('computing ... ' + string(j/Nt))
    % RTN rotation matrix
    rn_ACI = rn(:, j);
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    rn_ACAF = ACAF_ACI * rn_ACI;
    x = rn_ACAF(1); y = rn_ACAF(2); z = rn_ACAF(3); 
    lambda = atan2(y, x);  % in radians
    r = sqrt(x^2 + y^2);
    phi = atan2(z, r);

    ENU_ACAF = [
    -sin(lambda),             cos(lambda),              0;
    -sin(phi)*cos(lambda), -sin(phi)*sin(lambda),  cos(phi);
     cos(phi)*cos(lambda),  cos(phi)*sin(lambda),  sin(phi)
            ];

    % computed meas.
    [~, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm, ACAF_ACI'*ENU_ACAF');

     % compute Point mass approx error
    [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, ENU_ACAF');
    
    % % hc = Hc_RTN(9, 1:end);
    hx = Hp(9, 3);

    Z1 = Hc_RTN(9, 1:end);
    Z2 = [Hc_RTN(8, 1:end);Hc_RTN(7, 1:end)];

    hc = Z2;

    % compute information matrix
    Ax = Ax + hc' * inv(R) * hc;
    Mxc = Mxc + (hc' * inv(R) * hx);
end

% compute unertianty
Px = inv(Ax);
Pc = sigmar^2;
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';

sigma = sqrt(diag(Px));
sigmaC = sqrt(diag(Pxx));

sigma_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        sigma, Cnm.*0, Snm.*0); 
sigmaC_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        sigmaC, Cnm.*0, Snm.*0); 

[sigma2_Anly] = compute_gravDegreeVar(n_max, sigmam, GM, Re,...
    vecnorm(r0), Nt);
sigma_Anly = sqrt(sigma2_Anly);
sigma_Anly_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        sigma_Anly, Cnm.*0, Snm.*0); 

[sigma2_C_Anly, sigma2_Relerr] = compute_gravDegreeVar_PE(n_max, sigmam, sigmar, GM, Re,...
    vecnorm(r0), Nt, X);
sigma_C_Anly = sqrt(sigma2_C_Anly);
sigma_C_Anly_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        sigma_C_Anly, Cnm.*0, Snm.*0); 

% compute relative error
Rel_err =  sqrt((sigmaC.^2 -sigma.^2)./(sigma.^2));
sigma_err_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        Rel_err, Cnm.*0, Snm.*0);
sigma_err_Anly_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        sqrt(sigma2_Relerr), Cnm.*0, Snm.*0); 

figure()
semilogy(1:n_max, sigma_RMS, 'LineWidth', 2)
hold on;
semilogy(1:n_max, sigma_Anly_RMS, 'LineWidth', 2)
legend('numerical', 'analytical')
title('comparison between analytical and numerical grav. degree variance')

figure()
semilogy(1:n_max, sigma_RMS, 'LineWidth', 2)
hold all;
semilogy(1:n_max, sigmaC_RMS, 'LineWidth', 2)
semilogy(1:n_max, sigma_C_Anly_RMS, 'LineWidth', 2, ...
    'LineStyle', '--', 'Color', 'k')
legend('position free', 'numerical \sigma_r = 1', 'analytical')
title('comparison between position free and consider cov.')

figure()
semilogy(1:n_max, sigma_err_RMS, 'LineWidth', 2)
hold on;
semilogy(1:n_max, sigma_err_Anly_RMS, 'LineWidth', 2)
legend('numerical', 'analytical')
title('RMS Relative error coefficients variance')

% FUNCTIONS
function [sigma] = compute_DV(GM, Re, r, Nt, sigmaM, n)
    sigma = zeros(1, n);
    for j = 2:n
        %  %Z = (j+1) * (j+2);
        Z = (j+2) * sqrt(j * (j+1));
        sigma(j) = r^(3+j) / (Re^j) * 1 / GM * sigmaM / sqrt(Nt) * 1/ Z;
    end
end

function [sigma2] = compute_gravDegreeVar(n_max, sigmaM, GM, Re, r, Nt)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigma2 = ones(1, Ncs);
    alpha = 4*pi;
    n = 0;
    Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
    sigma2(1) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha);

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        sigma2(j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha);
        
        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        sigma2(j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha);

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [sigma2, sigma2_RelErr] = compute_gravDegreeVar_PE(n_max, sigmaM, sigmaR, GM, Re, r, Nt, X)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigma2 = ones(1, Ncs); sigma2_RelErr = sigma2;
    alpha = 4*pi;
    n = 0;
    Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
    An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
    sigma2(1) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
        (sigmaR * (An/Bn))^2 * X(1)^2;
    
    sigma2_RelErr(1) = Nt * (sigmaR/sigmaM)^2 * X(1)^2 * An^2;
    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        sigma2(j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
            (sigmaR * (An/Bn))^2 * X(j)^2;
        sigma2_RelErr(j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * An^2;
        
        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        sigma2(j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha);
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        sigma2(j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
            (sigmaR * (An/Bn))^2 * X(j)^2;
        sigma2_RelErr(j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * An^2;

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [sigmaK] = compute_Kaula(n_max, K)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigmaK = ones(1, Ncs);

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        sigmaK(j) = K/(n^2);
        
        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        % compute degree variance
         sigmaK(j) = K/(n^2);

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end