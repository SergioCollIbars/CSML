clear;
clc;
close all;
%%                             POST PROCESS
% ----------------------------------------------------------------------- %
%   Author: Sergio Coll Ibars
%
%   Date: 04/11/2022
%
%   Description: Post process program. Reads files from the data folder and
%   plots the results.
% ------------------------------------------------------------------------%

% Imports
addpath('../data_files/');
addpath('../modules/orbital_module/functions/');

% config plots
set(0,'defaultAxesFontSize',16);

% Input parameters
plotOrbit = 1;                                       % plot orbit boolean
plotAttitude = 0;                                    % plot Attitude boolean
plotAcc = 1;                                         % plot acc boolean
plotEstim = 1;                                       % plot estimation

test = false;                                        % show test plots

body = "Bennu";                                      % body name
path = "HARMCOEFS_BENNU_OSIRIS_0.txt";               % SH path
normalized = 0;                                      % normalized
d = 1;                                               % sensor distance [m]
G_constant = 6.6743E-11;                             % grav contant [m^3/Kg^2s^2]

if(body == "Bennu")
    rho = 1250;
    GM = rho * 4/3 * pi * (0.2449400088E3)^3 * G_constant;
    if(normalized == 0)
        GM = 5.2;
    end
elseif(body =="Eros")
    rho = 2670;
    GM = rho * 4/3 * pi * (8.507322E3)^3 * G_constant;
end

% start plots
if(plotAcc == true)
    % file name
    name = "accData.txt";

    % read data
    T = readtable(name);
        
    % extract data
    t = T.TIME';
    Nt = length(t);

    % create matrix
    ad = zeros(3, 3, length(t));
    Ub = zeros(3, length(t));
    
    % extract data
    for j = 1:Nt
        ad(1, 1, j) = T.ad_xx(j);
        ad(1, 2, j) = T.ad_xy(j);
        ad(1, 3, j) = T.ad_xz(j);

        ad(2, 1, j) = T.ad_yx(j);
        ad(2, 2, j) = T.ad_yy(j);
        ad(2, 3, j) = T.ad_yz(j);

        ad(3, 1, j) = T.ad_zx(j);
        ad(3, 2, j) = T.ad_zy(j);
        ad(3, 3, j) = T.ad_zz(j);
    end

    Ub(1, :) = T.Ub_x';
    Ub(2, :) = T.Ub_y';
    Ub(3, :) = T.Ub_z';

     % plot results
     AccP(t, ad, Ub, d);
end

if(plotOrbit == true)
    % file name
    name = "orbitData.txt";
    name2 = "N2SAR.txt";
    name3 = "N2ACAF.txt";

    % read data
    T = readtable(name);
    T2 = readtable(name2);
    T3 = readtable(name3);
        
    % extract data
    t = T.TIME';
    
    ri = zeros(3, length(t));
    ri(1, :) = T.ri_x;
    ri(2, :) = T.ri_y;
    ri(3, :) = T.ri_z;

    vi = zeros(3, length(t));
    vi(1, :) = T.vi_x;
    vi(2, :) = T.vi_y;
    vi(3, :) = T.vi_z;

    rb = zeros(3, length(t));
    rb(1, :) = T.rb_x;
    rb(2, :) = T.rb_y;
    rb(3, :) = T.rb_z;

    rSAR = zeros(3, length(t));
    vSAR = zeros(3, length(t));
    rACAF = zeros(3, length(t));
    vACAF = zeros(3, length(t));

    lat = zeros(1, length(t));
    lon = zeros(1, length(t));

    alphaSAR = zeros(8, length(t));

    for k = 1:length(t)
        up = k*3;
        down = up-2;
        SAR_N  = [T2.SAR_N_1(down:up), T2.SAR_N_2(down:up),...
            T2.SAR_N_3(down:up)];

        ACAF_N = [T3.ACAF_N_1(down:up), T3.ACAF_N_2(down:up), ...
            T3.ACAF_N_3(down:up) ];

        rSAR(:, k) = SAR_N * ri(:, k);
        vSAR(:, k) = SAR_N * vi(:, k);

        rACAF(:, k) = ACAF_N * ri(:, k);
        vACAF(:, k) = ACAF_N * vi(:, k);

        % compute longituda and latitude
        racaf = rACAF(:, k) / vecnorm(rACAF(:, k));
        lat(k) = atan2(racaf(3), sqrt(racaf(1)^2 + racaf(2)^2));
        lon(k) = atan2(racaf(2), racaf(1));

        % orbital elements in SAR frame
        [alphaSAR(:, k)] = orbitalElem(rSAR(:, k), vSAR(:, k), GM);
    end
    
    % orbital elements in N frame
    alpha = zeros(8, length(t));
    alpha(1, :) = T.e;
    alpha(2, :) = T.h;
    alpha(3, :) = T.a;
    alpha(5, :) = T.f;
    alpha(6, :) = T.i;
    alpha(7, :) = T.Omega;
    alpha(8, :) = T.omega;

    % plot resutls
    OrbitP(t, ri, rb, rACAF, lon, lat, body);

    if test == true
        % compute mean values
        alpha_m = zeros(7, length(t));

        % Orbital values
        J = 1.08E-3;
        R0 = 6563E3;
        mu = 3.9833e+14;

        % get initial values
        e0 = mean(alpha(1, :));
        a0 = mean(alpha(3, :));
        n0 = sqrt(mu / (a0^3));
        sigma0 = 0;
        i0 = mean(alpha(6, :));
        omega0 = alpha(8, 1);
        Omega0 = alpha(7, 1);

        for k = 1:length(t)
            % current time
            T = t(k);
            % compute mean values. J2 perturbation
            alpha_m(:, k) = J2_orbitalElem(R0, J, n0, a0, e0, i0, omega0, Omega0, sigma0, T);
        end
        
        % plot orbital elements comparison
        OrbitalElemP_comp(t, alpha, alpha_m);
    else
        % plot orbital elements. SAR frame
        OrbitalElemP(t, alphaSAR);
    end
    
end

if(plotAttitude == true)
    % file name
    name = "attitudeData.txt";

    % read data
    T = readtable(name);
        
    % extract data
    t = T.TIME';

    omega = zeros(3, length(t));
    Omega = zeros(3, length(t));

    theta = zeros(3, length(t));
    thetaDot = zeros(3, length(t));
    thetaDdot = zeros(3, length(t));

    omega(1, :) = T.omega_x;
    omega(2, :) = T.omega_y;
    omega(3, :) = T.omega_z;

    Omega(1, :) = T.Omega_x;
    Omega(2, :) = T.Omega_y;
    Omega(3, :) = T.Omega_z;
    
    theta(1, :) = T.theta1;
    theta(2, :) = T.theta2;
    theta(3, :) = T.theta3;

    thetaDot(1, :) = T.theta1Dot;
    thetaDot(2, :) = T.theta2Dot;
    thetaDot(3, :) = T.theta3Dot;

    thetaDdot(1, :) = T.theta1Ddot;
    thetaDdot(2, :) = T.theta2Ddot;
    thetaDdot(3, :) = T.theta3Ddot;

    % plot results
    AttitudeP(t, omega, Omega, theta, thetaDot, thetaDdot);

end

if(plotEstim == true)
    % file name
    name = "estimData_2.txt";
    name2 = "estimData_1.txt"; % postfit
    name3 = "estimData_3.txt"; % prefit
    name4 = "estimData_4.txt"; % covariance matrix
    name5 = "estimData_5.txt"; % covariance diag over time
    name6 = "estimData_6.txt"; % estate over time

    % read data
    T = readtable(name);
    T2 = readtable(name2);
    T3 = readtable(name3);
    T4 = readtable(name4);
    T5 = readtable(name5);
    T6 = readtable(name6);

    % obtain coefficients
    Xt = table2array(T6);
    X = T.X;
    TIME = T2.TIME;
    Nt = length(TIME);

    % obtain sigma value and coavariance value
    sigma = T.sigma;
    P = reshape(T4.P, [length(X), length(X)]);
    
    % obtain covariance diag over time
    sigma2_t = (table2array(T5)).^2;

    % obtain error
    err = T.CS_err;
    
    % coefficient number
    n_max = T.n_max(1);
    Nc = -2;
    for j = 1:n_max+1
        Nc = Nc + j;
    end
    Ns = Nc - n_max;
    
    % coefficient legend
    [coeffC, coeffS] = compute_coefLegend(n_max);

    % true coeff values
    [Cnm_Bennu, Snm_Bennu, R_Bennu] = readCoeff(path);
    [Xtrue] = mat2list(Cnm_Bennu, Snm_Bennu, Nc, Ns);

    Xtrue(1) = GM;

    Ncmax = Nc;
    Nsmax = Ns;

    % compute RSS truth values
    [RSS_truth] = computeRSS_coefErr(n_max, Nc, Ns, ...
        zeros(Nc+Ns, 1), Cnm_Bennu, Snm_Bennu, GM);

    [RSS_sigma] = computeRSS_coefErr(n_max, Nc, Ns, ...
       sigma, zeros(Nc+Ns), zeros(Nc+Ns), 0);

    % compute RMS coefficients
    Sc = sigma(1:Nc);                                       % sigma Cnm
    Ss = sigma(Nc+1:Ns+Nc);                                 % sigma Snm
    var = zeros(1, n_max);
    varZ = zeros(1, n_max);
    varST = zeros(1, n_max);

    % initialize values
    var(1) = Sc(1)^2;                                       % total variance
    varZ(1) = Sc(1)^2;                                      % zonal variance

    
    n = 2;
    m = 0;
    j = 2;
    sum = 0;
    while n <= n_max
        sum = sum + Sc(j)^2;
        if(m == 0)
            varZ(n) = Sc(j)^2;
        end

        if(m < n)
            m = m + 1;
        else
            var(n) = sum;
            varST(n) = sum - varZ(n);

            n = n + 1;
            m = 0;
            sum = 0;
        end
        j = j + 1;
    end

    n = 2;
    m = 1;
    j = 1;
    sum = 0;
    while n <= n_max
        sum = sum + Ss(j)^2;

        if(m < n)
            m = m + 1;
        else
            var(n) = var(n) + sum;
            varST(n) = varST(n) + sum;
            
            n = n + 1;
            m = 1;
            sum = 0;
        end
        j = j + 1;
    end
    
    % compute RSS
    STD_RSS = sqrt(var);
    disp("sigma RSS [n.d]")
    disp(STD_RSS);

    % compute RMS
    STD_RMS = zeros(1, n_max);
    STD_RMS_zonal = zeros(1, n_max);
    STD_RMS_ST = zeros(1, n_max);
    STD_RSS_zonal = zeros(1, n_max);
    STD_RSS_ST = zeros(1, n_max);

    STD_RMS(1) = sqrt(var(1));
    STD_RMS_zonal(1) = sqrt(varZ(1));
    STD_RMS_ST(1) = sqrt(varST(1));
    STD_RSS_zonal(1) = sqrt(varZ(1));
    STD_RSS_ST(1) = sqrt(varST(1));
    for n = 2:n_max
        STD_RMS(n) = sqrt(var(n) / (2*n + 1));
        STD_RMS_zonal(n) = sqrt(varZ(n));
        STD_RMS_ST(n) = sqrt(varST(n) / n);
        STD_RSS_zonal(n) = sqrt(varZ(n));
        STD_RSS_ST(n) = sqrt(varST(n));
    end
    disp("sigma RMS [n.d]")
    disp(STD_RMS);
  
    % plot Kaula bounding
    plotKaula(n_max, Cnm_Bennu, Snm_Bennu, STD_RMS_zonal,...
        STD_RMS_ST, normalized);

    % plot prefit && postfit
    plotPrefPosf(TIME, T2, T3);

    % compute RMS value postfit
    r = table2array(T2);
    for j = 1:6
        val = r(:, j+1);
        RMS = rms(val);
        disp("RMS dP_" + string(j) + " = " + string(RMS));
    end

    % compute correlation
    corr = zeros(Nc+Ns, Nc + Ns);
    for i = 1:Ns+Nc
        for j = 1:Ns+Nc
            corr(i, j) = abs(P(i, j) / (sqrt(P(i, i))*sqrt(P(j, j))));
        end
    end
    
    if(Nc ~= 1)
        figure();
        h1 = heatmap(coeffC,coeffC,corr(1:Nc, 1:Nc), 'Colormap',summer);
        title('Coefficient correlation, C_{nm}');

        figure();
        h2 = heatmap(coeffS,coeffS,corr(Nc+1:Ns+Nc, Nc+1:Ns+Nc), ...
            'Colormap',summer);
        title('Coefficient correlation, S_{nm}');
    end

       % Compute normalized estimation error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = sqrt((X(1) - GM)^2);
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefErr(n) = CoefErr(n)  + (X(j) - Cnm_Bennu(n + 1, m+1))^2;
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
        CoefErr(n) = CoefErr(n)  + (X(j) - Snm_Bennu(n + 1, m+1))^2;
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

    CoefTruth = zeros(1, n_max);
    CoefTruth(1) = sqrt((GM)^2);
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefTruth(n) = CoefTruth(n)  + (Cnm_Bennu(n + 1, m+1))^2;
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
        CoefTruth(n) = CoefTruth(n)  + (Snm_Bennu(n + 1, m+1))^2;
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
    for n = 2:n_max
        CoefTruth(n) = sqrt(CoefTruth(n) / (2*n + 1));
    end
    xAxis = linspace(0, n_max, n_max);
    xAxis(1) = 0;
    figure();
    semilogy(xAxis, CoefErr, 'r--sq', 'LineWidth', 1.5)
    hold on;
    semilogy(xAxis, CoefTruth, 'k-', 'LineWidth', 1.5)
    grid on;
    xlabel('Harmonic degree');
    ylabel('Coefficient error');
    title('RMS Coefficient error');
    legend('RMS error', 'RMS truth');

    figure();
    semilogy(xAxis, CoefTruth, 'k-', 'LineWidth', 1.5)
    grid on;
    hold on;
    semilogy(xAxis, STD_RMS, 'g--sq', 'LineWidth', 1.5)
    xlabel('Harmonic degree');
    ylabel('Coefficient error');
    title('RMS SH uncertainty, \sigma');
    legend('RMS truth', 'RMS \sigma');

    % plot statistics
    plotEstimation(Nc, Ns, Ncmax, X, Xtrue, err, ...
                    sigma, coeffC, coeffS)
    
    % plot gravity field map error
    plotMapAccError(X, n_max, Nc, Ns, Cnm_Bennu, Snm_Bennu, ...
        GM, R_Bennu, body)


    % plot covariance diag over time
    plotCovHist(sigma2_t, TIME, Nc, Ns)

    % plot extra parameters
    if(length(X) > (Nc + Ns))
        Ne = length(X) - (Nc + Ns);
        init = Nc + Ns + 1;
        txt = string(linspace(1, Ne, Ne));
        X0e = Xt(init:end, 1);

        % plot state over time
        figure()
        semilogy(TIME./86400, Xt(init:end, :), 'LineWidth', 1.5)
        xlabel('TIME [days]')
        ylabel('Value [-]')
        title('estate for extra parameters over time')
        legend(txt)

        % plot sigma
        figure()
        semilogy(TIME./86400, sqrt(sigma2_t(init:end, :)), 'LineWidth', 1.5)
        xlabel('TIME [days]')
        ylabel('\sigma value [-]')
        title('STD for extra parameters over time')
        legend(txt)

        % error in the estimation
        BT = [-282,8.88,16.42,-636,-19500,-113] * 1E-9 /2;
        B = [X0e(1), X0e(3), X0e(5), X0e(7), X0e(9), X0e(11)];
        DT = [0.0065, 0.0103, 0.0001, 0.195,0.9374,0.0033] * 1E-9/(2*86400);
        D = [X0e(2), X0e(4), X0e(6), X0e(8), X0e(10), X0e(12)];

        figure()
        xaxis = linspace(1, length(B), length(B));
        semilogy(xaxis, abs(BT), 'k-', 'LineWidth', 1.5)
        hold on;
        semilogy(xaxis, abs(B./1E10 - BT), 'b--sq', 'LineWidth', 1.5)
        title('Bias estimation error')
        legend('|Truth value|', 'Absolute error')
        xlabel('Bias')
        ylabel('error [1/s^2]')
        set(gca, 'XTick',1:length(xaxis), 'XTickLabel',...
            ["bxx", "bxy", "bxz", "byy", "byz", "bzz"])
        grid on;

        figure()
        xaxis = linspace(1, length(B), length(B));
        semilogy(xaxis, abs(DT), 'k-', 'LineWidth', 1.5)
        hold on;
        semilogy(xaxis, abs(D./1E10 - DT), 'b--sq', 'LineWidth', 1.5)
        title('Drift estimation error')
        legend('|Truth value|', 'Absolute error')
        xlabel('Drift')
        ylabel('error [1/s^3]')
        set(gca, 'XTick',1:length(xaxis), 'XTickLabel',...
            ["dxx", "dxy", "dxz", "dyy", "dyz", "dzz"])
        grid on;
    end
    
end


%%                          FUNCTIONS

function AccP(t, ad, Ub, d)
    % seconds 2 days
    t = t./86400;

    % acc per sensor
    figure();
    
    Ax = d;
    y1(:, :) = ad(1, 1, :)*Ax;
    y2(:, :) = ad(1, 2, :)*Ax;
    y3(:, :) = ad(1, 3, :)*Ax;

    subplot(4, 1, 1);
    plot(t, y1, t, y2, t, y3)
    xlabel('time [days]')
    ylabel('a_x [m / s^2]')
    grid on;
    legend('sensor x axis', 'sensor y axis', 'sensor z axis')

    y1(:, :) = ad(2, 1, :)*Ax;
    y2(:, :) = ad(2, 2, :)*Ax;
    y3(:, :) = ad(2, 3, :)*Ax;

    subplot(4, 1, 2);
    plot(t, y1, t, y2, t, y3)
    xlabel('time [days]')
    ylabel('a_y [m / s^2]')
    grid on;
    legend('sensor x axis', 'sensor y axis', 'sensor z axis')

    y1(:, :) = ad(3, 1, :)*Ax;
    y2(:, :) = ad(3, 2, :)*Ax;
    y3(:, :) = ad(3, 3, :)*Ax;

    subplot(4, 1, 3);
    plot(t, y1, t, y2, t, y3)
    xlabel('time [days]')
    ylabel('a_z [m / s^2]')
    grid on;
    legend('sensor x axis', 'sensor y axis', 'sensor z axis')

    y1(:, :) = sqrt((ad(1, 1, :)*Ax).^2 + (ad(2, 1, :)*Ax).^2 + (ad(3, 1, :)*Ax).^2);
    y2(:, :) = sqrt((ad(1, 2, :)*Ax).^2 + (ad(2, 2, :)*Ax).^2 + (ad(3, 2, :)*Ax).^2);
    y3(:, :) = sqrt((ad(3, 1, :)*Ax).^2 + (ad(2, 3, :)*Ax).^2 + (ad(3, 3, :)*Ax).^2);

    subplot(4, 1, 4);
    plot(t, y1, t, y2, t, y3)
    xlabel('time [days]')
    ylabel('|a| [m / s^2]')
    grid on;
    legend('sensor x axis', 'sensor y axis', 'sensor z axis')

    % tensor plot
    figure();
    
    lw = 1.5;
    colorval = "#D95319";   % phase B
    %colorval = "#EDB120";   % phase A
    subplot(3, 3, 1);
    scale = 1;
    y(:, :) = ad(1, 1, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{xx}')
    grid on;

    subplot(3, 3, 2);
    y(:, :) = ad(1, 2, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{xy}')
    grid on;

    subplot(3, 3, 3);
    y(:, :) = ad(1, 3, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{xz}')
    grid on;

    subplot(3, 3, 4);
    y(:, :) = ad(2, 1, :) * scale;
    plot(t, y,'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{yx}')
    grid on;

    subplot(3, 3, 5);
    y(:, :) = ad(2, 2, :) * scale;
    plot(t, y,'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{yy}')
    grid on;

    subplot(3, 3, 6);
    y(:, :) = ad(2, 3, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{yz}')
    grid on;

    subplot(3, 3, 7);
    y(:, :) = ad(3, 1, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{zx}')
    grid on;

    subplot(3, 3, 8);
    y(:, :) = ad(3, 2, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{zy}')
    grid on;

    subplot(3, 3, 9);
    y(:, :) = ad(3, 3, :) * scale;
    plot(t, y, 'LineWidth',lw, 'Color',colorval);
    xlabel('time [days]')
    ylabel('\Gamma_{zz}')
    grid on;
    sgtitle('Phase A. Gravity tensor components over time. Units [1/s^2]') 

    % Potential gradient plot. Body frame
    figure();

    subplot(3, 1, 1);
    plot(t, Ub(1, :));
    xlabel('time [days]');
    ylabel('Ub_x');
    grid on;
    legend('sensor x axis');
    title('Potential gradient diference. Body frame');

    subplot(3, 1, 2);
    plot(t, Ub(2, :, 1));
    xlabel('time [days]');
    ylabel('Ub_y');
    grid on;
    legend('sensor x axis');

    subplot(3, 1, 3);
    plot(t, Ub(3, :, 1));
    xlabel('time [days]');
    ylabel('Ub_z');
    grid on;
    legend('sensor x axis');

end

function OrbitP(t, ri, rb, rACAF, lon, lat, body)

% gray color
grayColor = [.7 .7 .7];

% plot 3D orbit
N = length(ri(1, :));
figure();
plot3(ri(1,:), ri(2,:), ri(3,:), 'LineWidth', 1.5);
hold on;
plot3(ri(1,:), ri(2,:),ones(1, N)*min(zlim), 'Color', grayColor);    % XY plane
plot3(ones(1, N)*min(xlim), ri(2,:), ri(3, :), 'Color',grayColor);   % ZY plane
plot3(ri(1, :), ones(1, N)*min(ylim), ri(3, :), 'Color', grayColor); % Zx plane
title('Orbit in the inertial frame {I, J, K}');
grid on;

% plot axis
mAxis = max(max(ri));
%axis([0 max(ri(1,:)) 0 max(ri(2,:)) 0 max(ri(3,:))])
axis([0 mAxis 0 mAxis 0 mAxis])
hold all;
quiver3(0,0,-max(0),0,0,max(zlim),'r','LineWidth',1)
quiver3(0,-max(0),0,0,max(ylim),0,'r','LineWidth',1)
quiver3(-max(0),0,0,max(xlim),0,0,'r','LineWidth',1)
text(0,0,max(zlim),'K','Color','r')
text(0,max(ylim),0,'J','Color','r')
text(max(xlim),0,0,'I','Color','r')

% create planet surface
[x,y,z] = sphere;

% Scale to desire radius.
if(body == "Earth")
    radius = 650000; % Earth radious [m]
    x = x * radius;
    y = y * radius;
    z = z * radius;

    % plot planet
    surf(x, y ,z);
    axis equal
elseif(body == "Bennu")
    scale =  450;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 

    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
elseif(body == "Eros")
    scale =  450;  % Eros object scale factor
    obj = readObj('Eros-poly.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 

    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end
view([145 20])

% plot scatter ACAF possition
figure();
scatter3(rACAF(1, :), rACAF(2, :), rACAF(3, :), 'LineWidth', 2);
title('Orbit in the ACAF frame');

% plot axis
mAxis = max(max(ri));
%axis([0 max(ri(1,:)) 0 max(ri(2,:)) 0 max(ri(3,:))])
axis([0 mAxis 0 mAxis 0 mAxis])
hold all;
quiver3(0,0,-max(0),0,0,max(zlim),'r','LineWidth',1)
quiver3(0,-max(0),0,0,max(ylim),0,'r','LineWidth',1)
quiver3(-max(0),0,0,max(xlim),0,0,'r','LineWidth',1)
text(0,0,max(zlim),'K','Color','r')
text(0,max(ylim),0,'J','Color','r')
text(max(xlim),0,0,'I','Color','r')

% create planet surface
[x,y,z] = sphere;

% Scale to desire radius.
if(body == "Earth")
    radius = 650000; % Earth radious [m]
    x = x * radius;
    y = y * radius;
    z = z * radius;

    % plot planet
    surf(x, y ,z);
    axis equal
elseif(body == "Bennu")
    scale =  450;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 

    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
elseif(body == "Eros")
    scale =  450;  % Eros object scale factor
    obj = readObj('Eros-poly.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 

    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end

% plot scatter surface
figure();
plot(rad2deg(lon), rad2deg(lat), '.');
grid on;
xlabel("Longitude, \lambda [deg]");
ylabel("Latitude, \Phi [deg]");
title('Planet sampling. ACAF frame');

% plot body position
figure();

subplot(4, 1, 1);
plot(t, rb(1, :));
xlabel('TIME [s]');
ylabel('r_b x component');
grid on;
title('Body frame position vector (r_b) ');

subplot(4, 1, 2);
plot(t, rb(2, :));
xlabel('TIME [s]');
ylabel('r_b y component');
grid on;

subplot(4, 1, 3);
plot(t, rb(3, :));
xlabel('TIME [s]');
ylabel('r_b z component');
grid on;

subplot(4, 1, 4);
plot(t, vecnorm(rb));
xlabel('TIME [s]');
ylabel('r_b magnitude');
grid on;

end

function plotEstimation(Nc, Ns, Ncmax, X, Xtrue, err, ...
                    sigma, coeffC, coeffS)
    % string cell length
    NcCoeff = length(coeffC);
    NsCoeff = length(coeffS);
    
    lw = 2;
    colorVal = "#D95319";   % phase B";

    % non dimension STD
    sigmaCU = sigma(1:Nc)./Xtrue(1:Nc);
    sigmaSU = sigma(Nc+1:Ns+Nc)./Xtrue(Ncmax+1:Ncmax+Ns);
    
    figure();
    
    errorbar(linspace(1,Nc-1, Nc-1), X(2:Nc), 3*sigma(2:Nc), ...
        'LineWidth', lw, 'color', colorVal)
    grid on;
    title('C_{nm} coefficient estimation error bar [n.d]')
    ylabel('C_{nm} +-3 \sigma')
    xlabel('nm order ');
    set(gca, 'XTick',1:NcCoeff-1, 'XTickLabel',coeffC(2:end))
    
    figure();
    
    errorbar(linspace(1,Ns, Ns), X(Nc+1:Nc+Ns), 3*sigma(Nc+1:Ns+Nc), ...
        'LineWidth', lw, 'color', colorVal)
    grid on;
    title('S_{nm} coefficient estimation error bar [n.d]')
    ylabel('S_{nm} +-3 \sigma')
    xlabel('nm order');
    set(gca, 'XTick',1:NsCoeff, 'XTickLabel',coeffS)

end

function plotMapAccError(X, n_max, Nxc, Nxs, Ct, St, GM, Re, body)
    % orbit parameters

    Np = 300;
    R = 1000;                           % Bennu SC radius

    % longitude meshgrid
    lat = linspace(-pi/2, pi/2, Np);
    lon = linspace(-pi, pi, Np);
    [Px, Py] = meshgrid(lon, lat);
    
    % compute Cmat and Smat at nominal
    C_mat = zeros(n_max + 1, n_max + 1);
    S_mat = zeros(n_max + 1, n_max + 1);

    C_mat(1, 1) = 1;
    
    n = 2;
    m = 0;
    for j = 2:Nxc
        N = n + 1;
        M = m + 1;
        C_mat(N, M) = X(j);
        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n +1;
        end
    end

    n = 2;
    m = 0;
    for j = Nxc + 1:Nxs + Nxc
        N = n + 1;
        M = m + 2;
        S_mat(N, M) = X(j);
        if(m < n - 1)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end
    
    % acc vector 
    acc_err = zeros(Np, Np);
    acc_t = acc_err;
    acc_e = acc_err;
    for j = 1:Np
        for i =1:Np
            % get current position. ACAF frame
            x = R * cos(Py(i, j)) * cos(Px(i, j));
            y = R * cos(Py(i, j)) * sin(Px(i, j));
            z = R * sin(Py(i, j));

            r_ACAF = [x;y;z];

            % compute estimated acceleration
            [~, dU2, ~] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                r_ACAF, Re, X(1), 0);
            % compute true acceleration
            [~, dU1, ~] = potentialGradient_nm(Ct, St, n_max, ...
                                                r_ACAF, Re, GM, 0);
            err = vecnorm(dU2 - dU1) / vecnorm(dU1) * 100;
            acc_t(i, j) = vecnorm(dU1);
            acc_e(i, j) = vecnorm(dU2);
            acc_err(i, j) = err;
        end
    end
    
    % plot 
    figure()
    h = pcolor(rad2deg(Px), rad2deg(Py), acc_t);
    hold on;
    contour(rad2deg(Px), rad2deg(Py), acc_t, 'k','LineWidth',1);
    colormap(jet);
    set(h, 'EdgeColor', 'none');
    xlabel('LONGITUDE [deg]');
    ylabel('LATITUDE [deg]')
    title('True acceleration value [m/s^2]')
    view(2);
    colorbar();

    
     figure()
    h = pcolor(rad2deg(Px), rad2deg(Py), acc_e);
    hold on;
    contour(rad2deg(Px), rad2deg(Py), acc_e, 'k','LineWidth',1);
    colormap(jet);
    set(h, 'EdgeColor', 'none');
    xlabel('LONGITUDE [deg]');
    ylabel('LATITUDE [deg]')
    title('Estimated acceleration value [m/s^2]')
    view(2);
    colorbar();
    

    figure()
    h = pcolor(rad2deg(Px), rad2deg(Py), acc_err);
    hold on;
    contour(rad2deg(Px), rad2deg(Py), acc_err, 'k','LineWidth',1);
    colormap(jet);
    set(h, 'EdgeColor', 'none');
    xlabel('LONGITUDE [deg]');
    ylabel('LATITUDE [deg]')
    title('Acceleration relative error')
    view(2);
    colorbar();
    
end

function plotPrefPosf(TIME, T2, T3)
    figure();
    r = table2array(T3);
    t = TIME./86400;
    for j = 1:6
        subplot(2, 3, j)
        plot(t, detrend(r(:, j+1)), 'r.')
        xlabel('TIME [days]')
        ylabel('dp_' + string(j));
        grid on;
    end
    sgtitle('Prefit values over time')

    figure();
    r = table2array(T2);
    t = TIME./86400;
    z = ["xx", "xy", "xz", "yy", "yz", "zz"];
    for j = 1:6
        subplot(2, 3, j)
        plot(t, r(:, j+1), 'b.')
        xlabel('TIME [days]')
        ylabel('P ' + z(j));
        grid on;
        [h, p] = chi2gof(r(:, j+1));
        disp('Goodness of fit ' + z(j) + '= ' + string(h) + ', ' + string(p))
    end
    sgtitle('Postfit values over time')
end

function plotKaula(n_max, Cnm_Bennu, Snm_Bennu, varZ, varST, normalized)
    % Kaula normalized bounds
    KaulaN_z = 0.084;
    KaulaN_s = 0.026;
    alphaN_s = 2.01;
    alphaN_z = 2.08;

    Kaula = zeros(2, n_max);

    for n= 1:n_max
        if(normalized == 1)
            N_zonal = 1;
            N_sectorial = 1;
        else
            [N_zonal] = NormFactor(n, 0);
            [N_sectorial] = NormFactor(n, 1);
        end

        Kaula(1, n) = (KaulaN_z/n^alphaN_z) * N_zonal;
        Kaula(2, n) = (KaulaN_s/n^alphaN_s) * N_sectorial;
    end

    % plot harmonics bounding
    zBound = zeros(1, n_max);
    sBound = zeros(1, n_max);
    zBound(1) = 1;
    sBound(1) = 0.026;
    xx = linspace(1, n_max, n_max);
    val = 0;
    n = 2;
    m = 0;
    while n  <= n_max
        N = n + 1;
        M = m + 1;
        if m == 0
            zBound(n) = abs(Cnm_Bennu(N, M));
        else
            val = val + Cnm_Bennu(N, M)^2 + Snm_Bennu(N, M)^2;
        end
        
        if(m < n)
            m = m + 1;
        else
            sBound(n) = sqrt(val / (2*n));
            val = 0;
            m = 0;
            n = n + 1;
        end
    end
     
    figure();

    subplot(1, 2, 1);
    semilogy(xx, zBound, 'r sq', xx, Kaula(1, 1:n_max), 'r --', ...
        xx, 3*varZ, 'g-sq', 'LineWidth', 1.5);
    grid on;
    xlabel('Hamonic order');
    ylabel('Uncertanty');
    legend('Truth Zonals', 'K_{zonal}', 'Zonal estimation 3 \sigma');

    subplot(1, 2, 2);
    semilogy(xx, sBound, 'r sq', xx, Kaula(2, 1:n_max), 'r --', ...
        xx, 3*varST, 'g-sq', 'LineWidth', 1.5)
    grid on;
    xlabel('Hamonic order');
    ylabel('Uncertanty');
    legend('Truth Sectorial', 'K_{sectorial}', ...
        'Sectorial estimation 3 \sigma');
    sgtitle('Truth harmonic bounding with estimation uncertanty');

end

function plotCovHist(sigma2_t, TIME, Nc, Ns)
    figure();
    
    subplot(2, 1, 1)
    semilogy(TIME./86400, sigma2_t(1, :), 'DisplayName',"GM");
    hold all;
    legend('-DynamicLegend');
    n = 2;
    m = 0;
    for nn = 2:Nc
        if(m == 0)
            semilogy(TIME./86400, sigma2_t(nn, :), 'DisplayName', ...
                "C_{" + string(n) + "0}");
        end
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    grid on;
    xlabel("TIME [days]");
    ylabel("1 \sigma^2")
    title('Zonal Cnm STD over time')

    subplot(2, 1, 2)
    n = 2;
    m = 0;
    for nn = 2:Nc
        if(m ~= 0)
            semilogy(TIME./86400, sigma2_t(nn, :), 'DisplayName', ...
                "C_{" + string(n) + string(m) + "}");
            legend('-DynamicLegend');
            hold all;
        end
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    grid on;
    xlabel("TIME [days]");
    ylabel("1 \sigma^2")
    title('Sectoral Cnm STD over time')

    figure();
    n = 2;
    m = 0;
    for nn = Nc+1:Nc+Ns
        semilogy(TIME./86400, sigma2_t(nn, :), 'DisplayName', ...
            "S_{" + string(n) + string(m+1) + "}");
        legend('-DynamicLegend');
        hold all;
 
        if(m < n-1)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    grid on;
    xlabel("TIME [days]");
    ylabel("1 \sigma^2")
    title('Snm STD over time')
end

function [CoefErr] = computeRMS_coefErr(n_max, Nc, Ns, X, Cnm, Snm, GM)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = sqrt((X(1) - GM)^2);
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

function [CoefErr] = computeRSS_coefErr(n_max, Nc, Ns, X, Cnm, Snm, GM)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = (X(1) - GM)^2;
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

function AttitudeP(t, omega, Omega, theta, thetaDot, thetaDdot)
% plot angular velocity
figure();

subplot(3, 1, 1);
plot(t, omega(1,:));
xlabel('Time [s]')
ylabel('\omega_x [rad/s]');
title('Angular velocity [rad/s], body frame')
grid on;

subplot(3, 1, 2);
plot(t, omega(2,:));
xlabel('Time [s]')
ylabel('\omega_y [rad/s]');
grid on;

subplot(3, 1, 3);
plot(t, omega(3,:));
xlabel('Time [s]')
ylabel('\omega_z [rad/s]');
grid on;

% plot angular acceleration
figure();

subplot(3, 1, 1);
plot(t, Omega(1,:));
xlabel('Time [s]')
ylabel('\Omega_x [rad/s]');
title('Angular acceleration [rad/s^2], body frame')
grid on;

subplot(3, 1, 2);
plot(t, Omega(2,:));
xlabel('Time [s]')
ylabel('\Omega_y [rad/s]');
grid on;

subplot(3, 1, 3);
plot(t, Omega(3,:));
xlabel('Time [s]')
ylabel('\Omega_z [rad/s]');
grid on;

% plot Euler profile
figure();

subplot(3, 1, 1);
plot(t, theta(1, :), t, theta(2, :), t, theta(3, :));
title('\theta profile');
xlabel('Time [s]');
legend('\theta_1', '\theta_2', '\theta_3');
grid on;

subplot(3, 1, 2);
plot(t, thetaDot(1, :), t, thetaDot(2, :), t, thetaDot(3, :));
title('\theta dot profile');
xlabel('Time [s]');
legend('\theta_1 dot', '\theta_2 dot', '\theta_3 dot');
grid on;

subplot(3, 1, 3);
plot(t, thetaDdot(1, :), t, thetaDdot(2, :), t, thetaDdot(3, :));
title('\theta ddot profile');
xlabel('Time [s]');
legend('\theta_1 ddot', '\theta_2 ddot', '\theta_3 ddot');
grid on;

end

function OrbitalElemP(t, alpha)
    % plot orbital elements
    t = t./86400; % time days
    lw = 2; % line width
    colorval = "#D95319";   % phase B
    %colorval = "#EDB120";   % phase A
    % Angular momentum magnitude
    figure();

    Avg = mean(alpha(2, :));
    epsilon = 1.005 * Avg;

    plot(t, alpha(2, :), 'Color', colorval);
    title('Angular momentum magnitude, h');
    xlabel('Time [s]');
    ylabel('h');
    grid on;
    ylim([Avg-epsilon, Avg+epsilon]);

    % orbital elements
    figure();
    
    Max = max(alpha(1, :));
    Min = min(alpha(1, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 1);
    plot(t, alpha(1, :), 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('e');
    grid on;
    ylim([Min, Max]);

    Max = max(alpha(3, :));
    Min = min(alpha(3, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 2);
    plot(t, alpha(3, :) , 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('a');
    grid on;
    ylim([Min, Max]);

    Max = max(alpha(5, :));
    Min = min(alpha(5, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 3);
    plot(t, rad2deg(alpha(5, :)), 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('f [deg]');
    grid on;
%     ylim([Min, Max]);
% 
%     Max = max(alpha(6, :));
%     Min = min(alpha(6, :));
%     if(Max == Min) 
%         Max = Max + 1E-6;
%     end

    subplot(3, 2, 4);
    plot(t, rad2deg(alpha(6, :)), 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('i [deg]');
    grid on;
%     ylim([Min, Max]);
% 
%     Max = max(alpha(7, :));
%     Min = min(alpha(7, :));
%     if(Max == Min) 
%         Max = Max + 1E-6;
%     end

    subplot(3, 2, 5);
    plot(t, rad2deg(alpha(7, :)), 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('\Omega [deg]');
    grid on;
%     if(~isnan(Max))
%         ylim([Min, Max]);
%     end
% 
%     Max = max(alpha(8, :));
%     Min = min(alpha(8, :));
%     if(Max == Min) 
%         Max = Max + 1E-6;
%     end

    subplot(3, 2, 6);
    plot(t, rad2deg(alpha(8, :)), 'LineWidth', lw, 'Color', colorval);
    xlabel('Time [days]');
    ylabel('\omega [deg]');
    grid on;
%     if(~isnan(Max))
%         ylim([Min, Max]);
%     end
    

    sgtitle('Orbital elements phase A. SAR frame');
end

function OrbitalElemP_comp(t, alpha, alpha_m)
    % plot orbital elements

    % Angular momentum magnitude
    figure();

    Avg = mean(alpha(2, :));
    epsilon = 1.005 * Avg;

    plot(t, alpha(2, :));
    title('Angular momentum magnitude, h');
    xlabel('Time [s]');
    ylabel('h');
    grid on;
    ylim([Avg-epsilon, Avg+epsilon]);

    % orbital elements
    figure();
    
    Max = max(alpha(1, :));
    Min = min(alpha(1, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 1);
    plot(t, alpha(1, :), t, alpha_m(1, :), '--');
    xlabel('Time [s]');
    ylabel('e');
    grid on;
    legend('e', 'e_m');
    ylim([Min, Max]);

    Max = max(alpha(3, :));
    Min = min(alpha(3, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 2);
    plot(t, alpha(3, :), t, alpha_m(2, :), '--');
    xlabel('Time [s]');
    ylabel('a');
    legend('a', 'a_m');
    grid on;
    ylim([Min, Max]);

    Max = max(alpha(5, :));
    Min = min(alpha(5, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 3);
    plot(t, alpha(5, :));
    xlabel('Time [s]');
    ylabel('f');
    grid on;
    ylim([Min, Max]);

    Max = max(alpha(6, :));
    Min = min(alpha(6, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 4);
    plot(t, alpha(6, :), t, alpha_m(3, :), '--');
    xlabel('Time [s]');
    ylabel('i');
    legend('i', 'i_m');
    grid on;
    ylim([Min, Max]);

    Max = max(alpha(7, :));
    Min = min(alpha(7, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 5);
    plot(t, alpha(7, :), t, alpha_m(4, :), '--');
    xlabel('Time [s]');
    ylabel('\Omega');
    grid on;
    legend('\Omega', '\Omega_m');
    ylim([Min, Max]);

    Max = max(alpha(8, :));
    Min = min(alpha(8, :));
    if(Max == Min) 
        Max = Max + 1E-6;
    end

    subplot(3, 2, 6);
    plot(t, alpha(8, :), t, alpha_m(5, :), '--');
    xlabel('Time [s]');
    ylabel('\omega');
    legend('\omega', '\omega_m');
    grid on;
    ylim([Min, Max]);
    

    sgtitle('Orbital elements comparison');
end

function [alpha] = J2_orbitalElem(R0, J2, n0, a0, e0, i0, omega0, Omega0, sigma0, t)
    % ------------------------------------------------------------------- %
    %                 COMPUTE J2 ORBITAL ELEMENTS FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 22/10/2022

    % Description: given the orbital element compute mean orbital elements
    %              for the J2 first order perturbation

    % Input:
    %   i: orbit inclination
    %   J2: planet shape constant
    %   R0: planet radious
    %   rho: orbital parameter
    %   e: orbit eccentricity
    %   t: time vector

    % Output:
    %   alpha: mean orbital elements vector, order:
        %   e_m: mean orbit eccentricity
        %   a_m: mean semi-axis 
        %   i_m: mean inclination
        %   Omega_m: mean ascending node
        %   omega_m: mean argument of periapsis
        %   n_m: mean mean motion
        %   sigma_m: mean nt

    % ------------------------------------------------------------------- %

    % compute orbital element at time step t.
    a_m = a0;
    e_m = e0;
    i_m = i0;
    
    i = i0;        
    e = e0;        
    rho = a_m * (1 - e_m^2);

    n_m = n0 * (1 + 3/2 * J2 * R0^2 / (rho^2) * (1 - 3/2 * sin(i)^2) * (1 - e^2)^(0.5));
    omega_m = omega0 + 3/2 * J2 * R0^2 / (rho^2) * n_m * (2 - 5 /2 * sin(i)^2) * t;
    Omega_m = Omega0 - 3/2 * J2 * R0^2 / (rho^2) * n_m * cos(i) * t;
    
    sigma_m = 3/2 * J2 * R0^2 / (rho^2) * (1 - 3/2 * sin(i)^2) * (1 - e^2)^(0.5)* n0 * t + sigma0;

    % Orbital elements vector 
    alpha = [e_m, a_m, i_m, Omega_m, omega_m, n_m, sigma_m];

end

function [Cmat, Smat, R] = readCoeff(path)
            list = table2array(readtable(path));
            degree = list(1);
            X = list(4:end);
            R = list(2);

            % count number of coefficients
            Nc = 1;
            for k = 2:degree
                Nc = Nc + k + 1;
            end
            Ns = 0;
            for k = 2:degree
                Ns = Ns + k;
            end
            
            % define matrices
            Cmat = zeros(degree + 1, degree + 1);
            Smat = zeros(degree + 1, degree + 1);
            Cmat(1, 1) = 1;

            n = 2;
            m = 0;
            for j = 2:Nc
                N = n + 1;
                M = m + 1;
       
                Cmat(N, M) = X(j);
                if(m < n)
                    m = m + 1;
                else
                    m = 0;
                    n = n +1;
                end
            end
        
            n = 2;
            m = 0;
            for j = Nc + 1:Ns + Nc
                N = n + 1;
                M = m + 2;
        
                Smat(N, M) = X(j);
                if(m < n - 1)
                    m = m + 1;
                else
                    m = 0;
                    n = n + 1;
                end
           
            end
end

function [X] = mat2list(C_mat, S_mat, Nxc, Nxs)
        X = zeros(Nxc + Nxs, 1);
        X(1) = 1;
        
        n = 2;
        m = 0;
        for j = 2:Nxc
            N = n + 1;
            M = m + 1;
            X(j) = C_mat(N, M);
            if(m < n)
                m = m + 1;
            else
                m = 0;
                n = n +1;
            end
        end

        n = 2;
        m = 0;
        for j = Nxc + 1:Nxs + Nxc
            N = n + 1;
            M = m + 2;
            X(j) = S_mat(N, M);
            if(m < n - 1)
                m = m + 1;
            else
                m = 0;
                n = n + 1;
            end
        end
end

function [N] = NormFactor(n, m)
    % Description: given degree, n and order, m compute the normalice
    % factor
    if(m == 0)
        delta = 1;
    else
        delta = 0;
    end
    fac1 = factorial(n - m);
    fac2 = factorial(n + m);
    N = ((2 - delta)*(2*n + 1) * fac1 /fac2)^(0.5);
end

function [coeffC, coeffS] = compute_coefLegend(n_max)
    % Description: compute the coefficient legend values
    
    coeffC = {'C00'};
    coeffS = {};

    n = 2;
    m = 0;
    
    while n <= n_max
        tt = 'C' + string(n) + string(m);
        coeffC{end + 1} = tt;

        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    n = 2;
    m = 1;

    while n <= n_max
        tt = 'S' + string(n) + string(m);
        coeffS{end + 1} = tt;
    
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end
 