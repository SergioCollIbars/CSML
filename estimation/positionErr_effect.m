clear;
clc;
close all;

%%          POSITION ERROR EFFECT IN GRAVITY ESTIMATION
% Description: Study the influence of position error in each SH coefficient
% for the simualted values in the GG.

% INCLUDE PATHS
addpath('functions/');

% READ DATA
T = readtable('accData.txt');
T2 = readtable('orbitData.txt');
dataAcc = table2array(readtable('accData.txt'));
dataOrbit = table2array(readtable('orbitData.txt'));

TIME = dataAcc(:, 1);
Nt = length(TIME);
At = TIME(2) - TIME(1);
G = dataAcc(:, 2:10);
ri = dataOrbit(:, 2:4)';
vi = dataOrbit(:, 8:10)';

% SH truth values
Cnm_t = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.0391557863539988 -2.96928723209235e-06 0.00375640654748657 0 0 0 0;...
        0.0148427177700986 0.00167095097673949 3.80845003468165e-05 0.000371755938456641 0 0 0;...
        0.0307494000000000 0.000413625917950024 -0.000490123739988179 -6.43092753252371e-05 4.51227856599555e-05 0 0;...
        -0.00456599734888228 0.000441961635589938 8.84264903209225e-05 1.40396087936725e-05 1.07458402471302e-05 1.19886311472876e-06 0;...
        -0.00896736657720649 0.000905916675449317 9.05780703119059e-05 2.77025499573633e-05 -1.92680576794137e-06 9.97533032070693e-08 1.67034838314692e-07];

Snm_t = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 0 -2.54325906400954e-05 0 0 0 0;...
        0 0.000992000134408593 6.53000000000000e-05 -0.00100797120329237 0 0 0;
        0 0.000634013000392474 0.000108054642426876 0.000102400000000000 0.00291093983173820 0 0;...
        0 -1.75754943031483e-05 -7.41878397813858e-05 4.79413750994751e-06 -0.000503800000000000 0.000448812426298560 0;...
        0 -1.35941163743731e-06 -9.69889209840526e-06 -7.55736000569855e-06 -1.59676058718751e-06 -0.000599300000000000 -3.93397896234722e-05];

% BIAS AND DRIFT CORRECTION VALUES. For phase B
    % Obtained
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
        
% Required
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
sigma2_ii = 6.32E-13;
sigma2_ij = 6.32E-13;
R0 = diag([sigma2_ii, sigma2_ij, sigma2_ij,...
            sigma2_ii, sigma2_ij, sigma2_ii]);
Ri = inv(R0);

Nm = 60;
R = linspace(0, 100, Nm);
T = zeros(1, Nm);
N = zeros(1, Nm);
deltaR = [R; N ;T];
Cerr_R = zeros(Nm, 6);
for j = 1:Nm
    disp('Iteriation = ' + string(j) + "/" + string(Nm))
    % LS ESTIMATOR
    Ax = 0;
    Nx = 0;
    for k = 1:Nt
        % ACI 2 ACAF rotation matrix
        Wt = W0 + W * TIME(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % current position. ECI frame
        r_ACI = ri(:, k);
        [ACI_RTN] = RTN2ECI(ri(:, k), vi(:, k));
        r_RTN = ACI_RTN' * r_ACI;

        % include error in RTN frame
        r_RTN = r_RTN + deltaR(:, j);
        
        % back to inertial
        r_ACI = ACI_RTN * r_RTN;
        r_ACAF = ACAF_ACI * r_ACI;
    
        % rotation matrix ECEF 2 body @ current time
        up = 3*k;
        down = up - 2;
        B_ACI = eye(3, 3);
        B_ACAF = B_ACI * ACAF_ACI';
        
        % get data measurements
        Ydata = -[G(k,1); G(k, 2); G(k, 3); G(k, 5); ...
            G(k, 6); G(k, 9)];
    
        deltaY = Ydata + B + k * D * At;
        
        % compute vibility matrix for the SH coefficients
        hi_9t = potentialGradient_Cnm(n_max, r_ACAF, ...
            Re, GM, B_ACAF, 0);
        H = [hi_9t(1, :); hi_9t(4, :); hi_9t(7, :); ...
             hi_9t(5, :); hi_9t(8, :); hi_9t(9, :)];
    
        % Batch filter
        [Ax, Nx] = batchFilter(deltaY, H, Ri, Ax, Nx);
    end
    
    % solve Normal equation
    XNOT = Ax\Nx;
    XNOT(1) = XNOT(1) * GM;
    
    % compute coefficient error
    Cerr_R(j, :) = computeRMS_coefErr(n_max, 26, 20, XNOT, Cnm_t, Snm_t);
    clc
end

% RMS SH coefficient
[RMS_SH] = computeRMS_coefErr(n_max, 26, 20, zeros(46, 1), Cnm_t, Snm_t);
Cerr1 =  Cerr_R;
Cerr2 =  Cerr_R;
for j = 1:6
    ind1 = find(Cerr_R(:, j) > RMS_SH(j));
    ind2 = find(Cerr_R(:, j) <= RMS_SH(j));
    Cerr1(ind1, j) = NaN;
    Cerr2(ind2, j) = NaN;
end

%% PLOT
figure()
semilogy(vecnorm(deltaR), Cerr_R, 'LineWidth', 1.5)
title('Position error effect in the estimation. Phase B')
xlabel('\Delta R [m]')
ylabel('RMS SH error')
legend('n = 1', 'n = 2', 'n = 3', 'n = 4', 'n = 5', 'n = 6')

figure()
semilogy(vecnorm(deltaR), Cerr1, 'LineWidth', 1.5)
hold on;
semilogy(vecnorm(deltaR), Cerr2, 'LineWidth', 1.5, 'Color','r')
title('Position error effect in the estimation. Phase B')
xlabel('\Delta R [m]')
ylabel('RMS SH error')
legend('n = 1', 'n = 2', 'n = 3', 'n = 4', 'n = 5', 'n = 6')

%% FUNCTION
function [CoefErr] = computeRMS_coefErr(n_max, Nc, Ns, X, Cnm, Snm)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = abs(X(1) - 5.2);
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

