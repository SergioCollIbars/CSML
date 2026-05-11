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

% Select body
body = "Earth";         % options: Bennu / Earth
if(body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
elseif(body == "Eros")
    % Asteroid parameters.
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM =  459604.431484721;
    n_max  = 10;
    normalized = 1;
else
    % Earth parameters
    path = "HARMCOEFS_EARTH_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 3.986004418E14;
    n_max  = 10;
    normalized = 1;
end
asterParams = [GM, Re, n_max, normalized];
n_max_K     = 240; % kaula max degree

% % % SH harmonics
[Nc, Ns, ~] = count_num_coeff(n_max_K); 
if(body == "Bennu")
    K = 0.025;
    [X_K] = compute_Kaula(n_max_K, K);
    X_K_RMS = computeRMS_coeffErr(n_max_K, Nc, Ns, ...
            X_K, zeros(n_max_K+1, n_max_K+1), zeros(n_max_K+1, n_max_K+1));
    X_K_RMS(1) = nan;
    
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    [X] = mat2list(Cnm, Snm, Nc, Ns);
    RMS_vals = computeRMS_coeffErr(n_max, Nc, Ns, ...
            X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

    X_RMS = X_K_RMS.*nan;
    X_RMS(1:length(RMS_vals)) = RMS_vals;

    % orbit radius
    r      = 1E3;
elseif(body == "Eros")
    K = 0.1;
    [X_K] = compute_Kaula(n_max_K, K);
    X_K_RMS = computeRMS_coeffErr(n_max_K, Nc, Ns, ...
            X_K, zeros(n_max_K+1, n_max_K+1), zeros(n_max_K+1, n_max_K+1));
    X_K_RMS(1) = nan;

    [Nc, Ns, ~] = count_num_coeff(n_max);
    [X] = mat2list(Cnm, Snm, Nc, Ns);
    RMS_vals = computeRMS_coeffErr(n_max, Nc, Ns, ...
            X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

    X_RMS = X_K_RMS.*nan;
    X_RMS(1:length(RMS_vals)) = RMS_vals;

    % orbit radius
    r      = 32E3;
else
    % compute Kaula rule
    K = 10^-5;
    [X_K] = compute_Kaula(n_max_K, K);
    X_K_RMS = computeRMS_coeffErr(n_max_K, Nc, Ns, ...
            X_K, zeros(n_max_K+1, n_max_K+1), zeros(n_max_K+1, n_max_K+1));
    X_K_RMS(1) = nan;

    [Nc, Ns, ~] = count_num_coeff(n_max);
    [X] = mat2list(Cnm, Snm, Nc, Ns);
    RMS_vals = computeRMS_coeffErr(n_max, Nc, Ns, ...
            X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

    X_RMS = X_K_RMS.*nan;
    X_RMS(1:length(RMS_vals)) = RMS_vals;

    % orbit radius
    r      = 350E3 + Re; 
end

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 9;
f = 1/1;
Nt = rev*T * f;

% specify number of meas assuming f = 1Hz
Nt = 31 * 86400 / 5;

% scale factor 
S = 1E6;

% weight matrix
sigmaM  = 3.3E-15;
sigmaR  = 1 / S;

% select measurement set
set = 1;

% compute degree variance analytical
[sigma2] = compute_gravDegreeVar(n_max_K, sigmaM, GM/(S^3), Re/S, r/S, Nt, set);
[Sxc] = compute_sensitivity(n_max, sigmaR, r/S, X);
[sigma2_r] = compute_gravDegreeVar_PE(n_max, sigmaM, GM/(S^3), Re/S, r/S, Nt, X, set);

Sxc_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            Sxc, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));  
sigmar_RMS  = computeRMS_coeffErr(n_max, Nc, Ns, ...
            sqrt(sigma2_r.*(S^2)), zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1)); 

[Nc, Ns, Ncs] = count_num_coeff(n_max_K); 
sigma_RMS  = computeRMS_coeffErr(n_max_K, Nc, Ns, ...
            sqrt(sigma2), zeros(n_max_K+1, n_max_K+1), zeros(n_max_K+1, n_max_K+1));

% plots
figure()
semilogy(1:n_max_K, sigma_RMS, 'LineWidth', 2)
hold all;
semilogy(1:n_max_K, X_K_RMS, 'LineWidth', 1.2, 'Color', 'k', 'LineStyle','--')
semilogy(1:n_max_K, X_RMS,   'LineStyle','none', 'Marker','*', 'Color','k')
xlabel('degree')
title('error degree variance')
legend('degree variance', 'Kaula')
grid on;
xlim([1, n_max_K])
% % figure()
% % semilogy(1:n_max, sigmar_RMS, 'LineWidth', 2, 'color', 'b')
% % xlabel('degree')
% % title('Position error threshold')
% % grid on;
% % 
% % figure()
% % semilogy(1:n_max, Sxc_RMS./X_RMS*100, 'LineWidth', 2, 'Color','g')
% % xlabel('degree')
% % ylabel('[%]')
% % title('Sensitivity error')
% % grid on;

%% FUNCTIONS
function [sigma2] = compute_gravDegreeVar(n_max, sigmaM, GM, Re, r, Nt, set)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigma2 = ones(1, Ncs);
    
    % n = 0
    n = 0;
    Z = selectSet(0, set);
    Bn = GM/(r^(3+n)) * Re^n;
    sigma2(1) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2;

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = GM/(r^(3+n)) * Re^n;
        Z = selectSet(n, set);
        sigma2(j) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2;
        
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
        Bn = GM/(r^(3+n)) * Re^n;
        Z = selectSet(n, set);
        sigma2(j) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2;

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [sigma2] = compute_gravDegreeVar_PE(n_max, sigmaM, GM, Re, r, Nt, X, set)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    sigma2 = ones(1, Ncs);
    n = 0;
    Z = selectSet(0, set);
    if(set == 1), F = 0; else, F = 2; end
    if(set == 4), F = 3; end
    Bn = GM/(r^(3+n)) * Re^n;
    sigma2(1, 1) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2 * (r/(3+n))^2  * 1/(F-1) * (1/ X(1))^2;

    
    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        Bn = GM/(r^(3+n)) * Re^n;
        Z = selectSet(n, set);
        sigma2(1, j) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2 * (r/(3+n))^2  * 1/(F-1) * (1/ X(j))^2;

        
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
        Bn = GM/(r^(3+n)) * Re^n;
        Z = selectSet(n, set);
        sigma2(1, j) = (1/Bn)^2 * (sigmaM^2/Nt) * (1/Z)^2 * (r/(3+n))^2  * 1/(F-1) * (1/ X(j))^2;

        % increment degree
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [Sxc] = compute_sensitivity(n_max, sigmaR, r, X)
    % define output
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    Sxc = ones(1, Ncs);
    n = 0;
    F = (1+n);
    Sxc(1) = (sigmaR * F / r) * X(1); % Ur,r

    % fill matrix
    m  = 0;
    n = 2;
    for j =2:Nc
        % compute degree variance
        F = (1+n);
        Sxc(j) = (sigmaR * F / r) * X(j); % Ur,r
        
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
        F = (1+n);
        Sxc(j) = (sigmaR * F / r) * X(j); % Ur,r

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

function [Z] = selectSet(n, set)
    ZMat = [(n+1)*(n+2); ...
        -(n+2)*sqrt(n*(n+1)); ...
        sqrt((n-1)*n*(n+1)*(n+2))];
    if(set == 4)
        Z = sqrt(ZMat(1)^2 + ZMat(2)^2);
    else
        Z = ZMat(set);
    end
end