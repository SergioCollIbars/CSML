clc;
close all;
clear;
format long g;

%%          ATTITUDE ERROR H MATRIX SENSITIVITY
% Description: Given some initial conditions, compute the sensitivity of 
%   H matrix to and attitude error.

% IMPUTS
addpath('functions/')

% INPUTS
measModel = 2; % 1: 9 meas / 2: 6 meas / 3: 5 meas
perturbation = "rnd"; % rnd = random / bias
N_mc  = 100;          % Monte Carlo samples

% Harmonics values
Cnm = [1 0 0 0 0 0 0;...
                    0 0 0 0 0 0 0;...
                    -0.0391557863539988 -2.96928723209235e-06 0.00375640654748657 0 0 0 0;...
                    0.0148427177700986 0.00167095097673949 3.80845003468165e-05 0.000371755938456641 0 0 0;...
                    0.0307494000000000 0.000413625917950024 -0.000490123739988179 -6.43092753252371e-05 4.51227856599555e-05 0 0;...
                    -0.00456599734888228 0.000441961635589938 8.84264903209225e-05 1.40396087936725e-05 1.07458402471302e-05 1.19886311472876e-06 0;...
                    -0.00896736657720649 0.000905916675449317 9.05780703119059e-05 2.77025499573633e-05 -1.92680576794137e-06 9.97533032070693e-08 1.67034838314692e-07];

Snm = [0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            0 0 -2.54325906400954e-05 0 0 0 0;...
            0 0.000992000134408593 6.53000000000000e-05 -0.00100797120329237 0 0 0;
            0 0.000634013000392474 0.000108054642426876 0.000102400000000000 0.00291093983173820 0 0;...
            0 -1.75754943031483e-05 -7.41878397813858e-05 4.79413750994751e-06 -0.000503800000000000 0.000448812426298560 0;...
            0 -1.35941163743731e-06 -9.69889209840526e-06 -7.55736000569855e-06 -1.59676058718751e-06 -0.000599300000000000 -3.93397896234722e-05];

% Bennu parameters
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
T = 2 * pi / n;
n_max = 6;
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% number of points
t0 = 0;
t_max = round(777600);
Nt_max = round(t_max / 5);
t = linspace(t0, t_max, Nt_max);

% Nominal orbit
X0 = [3.69978448523722;
    -105.053702653589;
    994.459667937083;
    -0.0528661190622108;
    -0.0487913496705586;
    -0.00495758536239498];
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM, Re, Cnm, Snm, n_max, 0, W, W0, ...
    RA, DEC), t, X0, options);
state = state';

% Position error. RTN
deltaR_RTN = [1; 1; 1];     % [m]
mean_delta = [0; 0; -0.5];  % [m]
sigma_delta = [0; 0; 1];    % [m]

% nominal orbit. ACI
r_ACI = zeros(3, Nt_max);
r_ACAF = zeros(3, Nt_max);
rr_ACI = zeros(3, Nt_max);

% Batch values
Ax = 0;
Nx = 0;
Axp = 0;
Nxp = 0;
sigma = 1E-40;
R = eye(9, 9) * sigma^2;

% MC loop
disp('MC run:')
errXp = zeros(N_mc, n_max);
tic
for j = 1:N_mc
    disp(j)
    % time loop
    for k = 1:Nt_max
        rt_ACI = [state(1, k);state(2, k);state(3, k)];
        vt_ACI = [state(4, k);state(5, k);state(6, k)];
        r_ACI(:, k) = rt_ACI;
        
        % ACAF postion. Truth
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        rt_ACAF = ACAF_ACI * rt_ACI;
        r_ACAF(:, k) = rt_ACAF;
    
        % RTN 2 ACI
        [ACI_RTN] = RTN2ECI(rt_ACI, vt_ACI);
        RTN_ACI = ACI_RTN';
    
        % RTN 2 ACAF
        ACAF_RTN = ACAF_ACI * ACI_RTN;
    
        % ACAF position. Perturbed
        if(perturbation == "bias")
            rp_ACAF = rt_ACAF + (ACAF_RTN * deltaR_RTN);
        elseif(perturbation == "rnd")
            delta = sigma_delta.*rand(size(mean_delta)) + mean_delta;
            rp_ACAF = rt_ACAF + (ACAF_RTN * delta);
        end
        rr_ACI(:, k) = ACAF_ACI' * rp_ACAF;
    
        % compute truth meas. RTN frame
        [~, ~, ddU] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                               rt_ACAF, Re, GM, 0);
        ddU = ddU + sigma * randn(size(ddU));
        Yo = reshape(ACAF_RTN'* ddU * ACAF_RTN, [9, 1]);
    
        % compute perturbed meas. RTN frame
        [~, ~, ddU] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                               rp_ACAF, Re, GM, 0);
        ddU = ddU + sigma * randn(size(ddU));
        Yc = reshape(ACAF_RTN'* ddU * ACAF_RTN, [9, 1]);
    
        % H truth. RTN
        [Ht] = potentialGradient_Cnm(n_max, rt_ACAF, Re, GM, ACAF_RTN');
    
        % H perturbed. ENU
        [Hp] = potentialGradient_Cnm(n_max, rp_ACAF, Re, GM, ACAF_RTN');
    
        % new meas
        if(measModel == 2)
            Y = Yo;
            Y = [Y(1); Y(2); Y(3); Y(5); Y(6); Y(9)];
            Ht = [Ht(1,:); Ht(2, :); Ht(3, :); Ht(5,:); Ht(6,:); Ht(9, :)];
            Hp = [Hp(1,:); Hp(2, :); Hp(3, :); Hp(5,:); Hp(6,:); Hp(9, :)];
            R = eye(6, 6) * sigma^2;
        elseif(measModel == 3)
            Y = Yo - Yc;
            Y = [Y(1) - Y(5); Y(2); Y(3); Y(6); Y(9)];
            Ht = [Ht(1,:) - Ht(5,:); Ht(2, :); Ht(3, :); Ht(6,:); Ht(9, :)];
            Hp = [Hp(1,:) - Hp(5,:); Hp(2, :); Hp(3, :); Hp(6,:); Hp(9, :)];
            R = eye(5, 5) * sigma^2;
        end
    
        % error. ENU coords
        H_err = Ht - Hp;
    
        Ax = Ax + (Ht' * inv(R) * Ht);
        Nx = Nx + (Ht' * inv(R) * Y);
    
        Axp = Axp + (Hp' * inv(R) * Hp);
        Nxp = Nxp + (Hp' * inv(R) * Y);
    end

    % Solve error
    X = Ax\Nx;
    Xp = Axp\Nxp;
    X(1) = GM * X(1);
    Xp(1) = GM * Xp(1);
    
    Px = sqrt(inv(Ax));
    Pxp = sqrt(inv(Axp));
    
    % compute variance RMS
    if(n_max == 4)
        Nc = 13;
        Ns = 9;
    elseif(n_max == 6)
        Nc = 26;
        Ns = 20;
    end
    
    % compute relative error. RMS
    [RMS_C_S] = computeRMS_coefErr(n_max, Nc, Ns, zeros(Nc+Ns, 1), Cnm, Snm);
    [errX] = computeRMS_coefErr(n_max, Nc, Ns, X, Cnm, Snm);
    [errXp(j, :)] = computeRMS_coefErr(n_max, Nc, Ns, Xp, Cnm, Snm);
    err1 = abs(errX - errXp(j, :));
    
    % compute relative error. RSS
    [RSS_C_S] = computeRSS_coefErr(n_max, Nc, Ns, zeros(Nc+Ns, 1), Cnm, Snm);
    [RSS_errX] = computeRSS_coefErr(n_max, Nc, Ns, X, Cnm, Snm);
    [RSS_errXp] = computeRSS_coefErr(n_max, Nc, Ns, Xp, Cnm, Snm);
    err2 = abs(errX - errXp);

end
toc

% plot orbit
figure();
subplot(1, 2, 1)
plot3(r_ACAF(1,:), r_ACAF(2, :), r_ACAF(3, :), 'LineWidth', 1.5);
axis equal;
grid on;
title('ACAF frame')

subplot(1, 2, 2)
plot3(r_ACI(1,:), r_ACI(2, :), r_ACI(3, :), 'LineWidth', 1.5);
hold on;
plot3(rr_ACI(1,:), rr_ACI(2, :), rr_ACI(3, :), 'LineWidth', 1.5);
axis equal;
grid on;
legend('Truth', 'Perturbed');
title('ACI frame')

if(perturbation == "bias")
    figure()
    plot3(r_ACI(1,:), r_ACI(2, :), r_ACI(3, :), 'LineWidth', 1.5);
    hold on;
    plot3(rr_ACI(1,:), rr_ACI(2, :), rr_ACI(3, :), 'LineWidth', 1.5);
    axis equal;
    grid on;
    legend('Truth', 'Perturbed');
    title('ACI frame')
elseif(perturbation == "rnd")
    figure()
    plot3(r_ACI(1,:), r_ACI(2, :), r_ACI(3, :), 'LineWidth', 1.5);
    hold on;
    scatter3(rr_ACI(1,:), rr_ACI(2, :), rr_ACI(3, :), 'SizeData', 1.5);
    axis equal;
    grid on;
    legend('Truth', 'Perturbed');
    title('ACI frame')
end

nval = linspace(1, n_max, n_max);
if(perturbation == "bias")
    % covariance stats
    [sx] = computeRSS_variance(diag(Px), Nc, n_max);
    [sxp] = computeRSS_variance(diag(Pxp), Nc, n_max);
    disp('RMS variance:')
    disp(sx);
    disp('RMS variance % error:')
    disp((sx - sxp)./sx * 100)

    % plot error. RMS
    figure();
    subplot(2, 1, 1);
    semilogy(nval, errX, 'r--sq', 'LineWidth', 1.5);
    hold on;
    semilogy(nval, errXp, 'b--sq', 'LineWidth', 1.5);
    hold on;
    semilogy(nval, RMS_C_S, 'k-sq', 'LineWidth', 1.5);
    xlabel('Harmonic degree, n');
    ylabel('Coefficient error');
    title('RMS coefficient error. \sigma = 10^{-40} 1/s^2');
    grid on;
    legend('Nominal', 'Perturbed', 'Truth coefficient');
    
    subplot(2, 1, 2);
    bar(nval, err1, 'cyan')
    xlabel('Harmonic degree, n');
    ylabel('RMS error difference');
    title('Difference between perturbed and truth error');
    grid on;
    
    % plot error. RSS
    figure();
    subplot(2, 1, 1);
    semilogy(nval, RSS_errX, 'r--sq', 'LineWidth', 1.5);
    hold on;
    semilogy(nval, RSS_errXp, 'b--sq', 'LineWidth', 1.5);
    hold on;
    semilogy(nval, RSS_C_S, 'k-sq', 'LineWidth', 1.5);
    xlabel('Harmonic degree, n');
    ylabel('Coefficient error');
    title('RSS coefficient error. \sigma = 10^{-40} 1/s^2');
    grid on;
    legend('Nominal', 'Perturbed', 'Truth coefficient');
    
    subplot(2, 1, 2);
    bar(nval, err2, 'cyan')
    xlabel('Harmonic degree, n');
    ylabel('RSS error difference');
    title('Difference between perturbed and truth error');
    grid on;
elseif(perturbation == "rnd")
    % covariance stats
    disp('RMS error variance:')
    disp(var(errXp));
    disp('RMS error mean');
    disp(mean(errXp));

    % coefficient RMS error MC 
    figure();
    scatter(nval, errXp);
    set(gca,'yscale','log')
    hold on;
    semilogy(nval, RMS_C_S, 'k-sq', 'LineWidth', 1.5);
    xlabel('Harmonic degree, n');
    ylabel('Coefficient error');
    title("RMS coefficient error MC ="  + string(N_mc) + " samples. " + ...
        "\sigma =" + string(sigma) + " 1/s^2");
    grid on;
    legend('Perturbed', 'Truth coefficient');
end


%% FUNCTIONS
function [s] = computeRSS_variance(P, Nc, n_max)
    % Description: compute RSS (root sum square) of the variance
    s = zeros(1, n_max);
    s(1) = P(1)^2;

    n = 2;
    m = 0;
    j = 2;
    while n <= n_max
        s(n) = s(n) + P(j)^2;
        if(m < n)
            m = m + 1;
        else 
            n = n + 1;
            m = 0;
        end
        j = j + 1;
    end
    n = 2;
    m = 1;
    j = Nc+1;
    while n <= n_max
        s(n) = s(n) + P(j)^2;
        if(m < n)
            m = m + 1;
        else 
            n = n + 1;
            m = 1;
        end
        j = j + 1;
    end
end


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

function [CoefErr] = computeRSS_coefErr(n_max, Nc, Ns, X, Cnm, Snm)
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
    for n = 1:n_max
        CoefErr(n) = sqrt(CoefErr(n));
    end
end



