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

% % % Asteroid parameters.
% % path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % GM = 5.2;
% % n_max  = 6;
% % normalized = 1;

% Earth parameters
savedData = 0;                % use saved data. 1 = yes / 0 = no
path = "HARMCOEFS_EARTH_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
path = "SIGMACOEFS_EARTH_1.txt";
[sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
GM = 3.986004418E14;
n_max  = 200;
normalized = 1;

asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
% % [X] = mat2list(Cnm, Snm, Nc, Ns);
% % X_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
% %         X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

% compute Kaula rule
% % K = 0.026;
K = 10^-5;
[CS_K] = compute_Kaula(n_max, K);
CS_K_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        CS_K, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

% Initial conditions
% % r      = 0.35E3;
r      = 250E3 + Re; 

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 1*16;
f = 1/1;
Nt = rev*T * f;

% scale factor 
S = 1E6;

% weight matrix
sigmaM  = 1E-12;
sigmaR  = 2E-2 / S;

% compute degree variance analytical
[sigma2] = compute_gravDegreeVar(n_max, sigmaM, GM/(S^3), Re/S, r/S, Nt);
[sigma2_PE, sigma2_RelErr] = compute_gravDegreeVar_PE(n_max, sigmaM, sigmaR, GM/(S^3), Re/S, r/S, Nt, CS_K);
[Sxc] = compute_sensitivity(n_max, sigmaR, GM/(S^3), Re/S, r/S, ones(1, Ncs));
sigma_PE_RMS = zeros(2, n_max); sigma_Err_RMS = sigma_PE_RMS;
for j = 1:2
    val  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            sqrt(sigma2_PE(j, :)), zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
    sigma_PE_RMS(j, :) = val;
    val  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            sqrt(sigma2_RelErr(j, :)), zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
    sigma_Err_RMS(j, :) = val;
end
Sxc_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            Sxc, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 
sigma_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            sqrt(sigma2), zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 

% plots
figure()
semilogy(1:n_max, sigma_RMS, 'LineWidth', 2)
hold all;
semilogy(1:n_max, sigma_PE_RMS, 'LineWidth', 2)
semilogy(1:n_max, CS_K_RMS, 'LineWidth', 2, 'Color', 'k')
xlabel('order')
title('RMS order variance, \sigma_r = ' + string(sigmaR*S) + ' m')
legend('No pos error', 'U_{rr}', 'U_{r \phi}', 'Kaula')
grid on;

figure()
semilogy(1:n_max, sigma_Err_RMS.*100, 'LineWidth', 2)
xlabel('order')
title('RMS variance relative error, \sigma_r = ' + string(sigmaR*S) + ' m')
legend('U_{rr}', 'U_{r \phi}')
ylabel('[%]')
grid on;

figure()
semilogy(1:n_max, Sxc_RMS.*100, 'LineWidth', 2, 'Color','g')
xlabel('order')
ylabel('[%]')
title('Position error sensitivity')
grid on;

figure()
semilogy(1:n_max, sigma_RMS./X_RMS.*100, 'LineWidth', 2)
hold all;
semilogy(1:n_max, sigma_PE_RMS(1, :)./X_RMS.*100, 'LineWidth', 2)
xlabel('order')
title('NSR RMS value, \sigma_r = ' + string(sigmaR*S) + ' m')
legend('No pos error', 'U_{rr}')
grid on;
ylabel('[%]')

%% FUNCTIONS
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
    sigma2 = ones(2, Ncs); sigma2_RelErr = sigma2;
    alpha = 4*pi;
    n = 0;
    Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
    An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
    Bn2 = sqrt(n) * sqrt(n+1) * (n+2) *GM/(r^(3+n)) * Re^n;
    An2 = (3+n)*(2+n)*GM/(r^(4+n)) * Re^n;
    sigma2(1, 1) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
        (sigmaR * (An/Bn))^2 * X(1)^2; % Ur, r
    sigma2(2, 1) =(sigmaM/Bn2)^2 *(4*pi/Nt) * (1/ alpha) + ...
        (sigmaR * (An2/Bn2))^2 * X(1)^2; % Ur, phi
    
    sigma2_RelErr(1, 1) = Nt * (sigmaR/sigmaM)^2 * X(1)^2 * An^2;     % Ur, r
    sigma2_RelErr(2, 1) = Nt * (sigmaR/sigmaM)^2 * X(1)^2 * (An2)^2 * n * (n+1); % Ur, phi
    
    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        Bn2 = sqrt(n) * sqrt(n+1) * (n+2) *GM/(r^(3+n)) * Re^n;
        An2 = (3+n)*(2+n)*GM/(r^(4+n)) * Re^n;
        sigma2(1, j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
            (sigmaR * (An/Bn))^2 * X(j)^2; % Ur, r
        sigma2(2, j) = (sigmaM/(Bn2))^2 * (4*pi/Nt) * (1/ alpha) + ...
        (sigmaR * (An2/Bn2))^2 * X(j)^2;     % Ur, phi
        sigma2_RelErr(1, j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * An^2;     % Ur, r
        sigma2_RelErr(2, j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * (An2)^2 * n * (n+1); % Ur, phi
        
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
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        Bn2 = sqrt(n) * sqrt(n+1) * (n+2) *GM/(r^(3+n)) * Re^n;
        An2 = (3+n)*(2+n)*GM/(r^(4+n)) * Re^n;
        sigma2(1, j) = (sigmaM/Bn)^2 * (4*pi/Nt) * (1/ alpha) + ...
            (sigmaR * (An/Bn))^2 * X(j)^2;  % Ur, r
        sigma2(2, j) = (sigmaM/(Bn2))^2 * (4*pi/Nt) * (1/ alpha) + ...
        (sigmaR * (An2/Bn2))^2 * X(j)^2;      % Ur, phi
        sigma2_RelErr(1, j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * An^2; % Ur, r
        sigma2_RelErr(2, j) = Nt * (sigmaR/sigmaM)^2 * X(j)^2 * (An2)^2 * n * (n+1); % Ur, phi

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [Sxc] = compute_sensitivity(n_max, sigmaR, GM, Re, r, X)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    Sxc = ones(1, Ncs);
    n = 0;
    Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
    An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
    Sxc(1) = (sigmaR * (An/Bn)) * X(1); % Ur,r

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = (2+n)*(1+n)*GM/(r^(3+n)) * Re^n;
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        Sxc(j) = (sigmaR * (An/Bn)) * X(j);
        
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
        An = (3+n)*(2+n)*(1+n)*GM/(r^(4+n)) * Re^n;
        Sxc(j) = (sigmaR * (An/Bn)) * X(j);

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