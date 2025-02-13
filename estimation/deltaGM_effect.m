clear;
clc;
close all;

%%          POSITION EFFECT ERROR IN GRAVITY ESTIMATION
% Description: Study the influence of GM error in each SH coefficeint
% for the simualted values in the GG.

% INCLUDE PATHS
addpath('functions/');

% READ DATA
T = readtable('accData.txt');
T2 = readtable('orbitData.txt');
dataAcc = table2array(readtable('accData.txt'));
dataOrbit = table2array(readtable('orbitData.txt'));
path = "HARMCOEFS_BENNU_OSIRIS_0.txt";
[Cnm_t, Snm_t] = readCoeff(path);

TIME = dataAcc(:, 1);
N = length(TIME);
At = TIME(2) - TIME(1);
G = dataAcc(:, 2:10);
ri = dataOrbit(:, 2:4)';
vi = dataOrbit(:, 8:10)';


% BIAS AND DRIFT CORRECTION VALUES. For phase B
% %    % Actual
% % B = [-1.41000139962798e-07;...
% % 4.42692393515439e-09;...
% % 8.21383297049534e-09;...
% % -3.18001688021536e-07;...
% % -9.74999439415953e-06;...
% % -5.65000358162869e-08];
% % 
% % D = [ 3.78130353377026e-17;...
% % 8.11853407415446e-17;...
% % -6.23503707619306e-18;...
% % 1.13081945673253e-15;...
% % 5.41557860429089e-15;...
% % 1.96147289639821e-17];
        
% Ideal
B = [-1.40999980249042e-07;...
4.44012601351793e-09;...
8.20964327520085e-09;...
-3.18000392726029e-07;...
-9.75000015292926e-06;...
-5.64996665239974e-08];

D = [3.75113030820805e-17;...
5.94765447757001e-17;...
7.04643715521234e-19;...
1.12864253647863e-15;...
5.42507893184713e-15;...
1.90770566034683e-17];

% % % bias value. Truth
% % b = [-282,8.88,16.42;1848,-636,-19500;5,9699,-113] ...
% % * 1E-9/2;  % [1/s^2]
% % 
% % 
% % % drift value. Truth
% % d = [0.0065,0.0103, 0.0001; 3.123,0.195,0.9374;0.001, 2.98,0.0033] ...
% % * 1E-9/(2*86400); % [1/s^2/s]
% % 
% % B = [b(1,1);b(1,2);b(1,3);b(2,2);b(2,3);b(3,3)];
% % D = [d(1,1);d(1,2);d(1,3);d(2,2);d(2,3);d(3,3)];

% ASTEROID PARAMETERS. Bennu
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
n_max = 6;
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% NOISE VALUE
sigma2_ii = 8*6.32E-13;
sigma2_ij = 8*6.32E-13;
Ri = diag([sigma2_ii, sigma2_ij, sigma2_ij,...
            sigma2_ii, sigma2_ij, sigma2_ii].^2);

Nm = 10;
deltaGM = linspace(-1, 1, Nm);
Cerr_GM = zeros(Nm, 6);
for j = 1:Nm
    disp('Iteriation = ' + string(j) + "/" + string(Nm))
    % LS ESTIMATOR
    Ax = 0;
    Nx = 0;
    for k = 1:N
        % ACI 2 ACAF rotation matrix
        Wt = W0 + W * TIME(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % current position. ECI frame
        r_ACI = ri(:, k);
        r_ACAF = ACAF_ACI * r_ACI;
    
        % rotation matrix ECEF 2 body @ current time
        up = 3*k;
        down = up - 2;
        B_ACI = eye(3, 3);
        B_ACAF = B_ACI * ACAF_ACI';
        
        % get data measurements
        [~, ~, ddU] = potentialGradient_nm(Cnm_t, Snm_t, 0, ...
                                                r_ACAF, Re, GM + deltaGM(j), 0);
        ddU = B_ACAF * ddU * B_ACAF';
        y = [ddU(1,1);ddU(1,2);ddU(1,3);ddU(2,2);ddU(2,3);ddU(3,3)];

        Ydata = -[G(k,1); G(k, 2); G(k, 3); G(k, 5); ...
            G(k, 6); G(k, 9)];
    
        deltaY = Ydata + B + k * D * At - y;
        
        % compute vibility matrix for the SH coefficients
        [~, hi_9t] = potentialGradient_Cnm(n_max, r_ACAF, ...
            Re, GM + deltaGM(j), B_ACAF, 0);
        hi_9t(:, 1) = hi_9t(:, 1) / (GM + deltaGM(j)); 
        H = [hi_9t(1, :); hi_9t(4, :); hi_9t(7, :); ...
             hi_9t(5, :); hi_9t(8, :); hi_9t(9, :)];
    
        % Batch filter
        [Ax, Nx] = batchFilter(deltaY, H, Ri, Ax, Nx);
    end
    
    % solve Normal equation
    XNOT = Ax\Nx;
    XNOT(1) = XNOT(1) + (GM + deltaGM(j));
    
    % compute coefficient error
    Cerr_GM(j, :) = computeRMS_coefErr(n_max, 26, 20, XNOT, Cnm_t, Snm_t);
    clc
end

% RMS value coefficients
RMS_val = computeRMS_coefErr(n_max, 26, 20, XNOT.*0, Cnm_t, Snm_t);
RMS_val(1) = sqrt(RMS_val(1));

%% PLOT
figure()
plot(deltaGM./GM, Cerr_GM./RMS_val * 100, 'LineWidth', 1.5)
hold on;
%plot(ones(1, Nm) * 5.2, linspace(0, 0.16, Nm), 'r', 'LineWidth', 1.5)
title('GM error effect in the estimation. Phase B')
xlabel('\delta GM / GM [-]')
ylabel('% error')
legend('GM', 'n = 2', 'n = 3', 'n = 4', 'n = 5', 'n = 6')

%% FUNCTION
function [CoefErr] = computeRMS_coefErr(n_max, Nc, Ns, X, Cnm, Snm)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = (X(1) - 5.2)^2;
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefErr(n) = CoefErr(n)  + (X(j) - Cnm(n + 1, m+1))^2;
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
        CoefErr(n) = CoefErr(n)  + (X(j) - Snm(n + 1, m+1))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
    for n = 2:n_max
        CoefErr(n) = sqrt(CoefErr(n) / (2*n + 1));
    end
end

