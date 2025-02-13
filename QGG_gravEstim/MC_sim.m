clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

%%              MONTE CARLO SIMULATION
% Inputs
MC = 200;         % number of MC runs

% outputs
CoefErr_RMS = zeros(MC, 6);
CoefErr = zeros(MC, 46);

% run MC simulation
for mc = 1:MC
    disp("MC simulation: " + string(mc) + "/" + string(MC))
    % run main
    run("main.m")

    % read estimation results
    T = readtable("estimData_2.txt");
    X = T.X;
    err = T.CS_err;
    CoefErr(mc, :) = err(1:46);
    [CoefErr_RMS(mc, :)] = computeRMS_coefErr(6, 26, 20, err);
    clc;
end
% compute consider covariance
T4 = readtable("estimData_4.txt");
T6 = readtable("estimData_6.txt");
Px = reshape(T4.P, [length(X), length(X)]);
Pc = Px;

% compute 1 sigma bounding
sigma = T.sigma;
sigmac = sqrt(diag(Pc));
bound_up = sigma' + median(CoefErr); 
bound_down = - sigma' + median(CoefErr); 
boundC_up = sigmac' + median(CoefErr); 
boundC_down = - sigmac' + median(CoefErr); 

sigma_RMS = computeRMS_coefErr(6, 26, 20, sigma);
sigmac_RMS = computeRMS_coefErr(6, 26, 20, sigmac);
boundRMS_up = std(CoefErr_RMS) + mean(CoefErr_RMS); 
boundRMS_down = - (std(CoefErr_RMS) + mean(CoefErr_RMS)); 
boundCRMS_up = sigmac_RMS + median(CoefErr_RMS); 
boundCRMS_down = - sigmac_RMS + median(CoefErr_RMS); 

% plot RMS error vs truth
truth = [5.2, 0.0175913977026888, 0.00568219239266755, 0.0102690380080387, ...
    0.00138487808268537, 0.00249992481290561];
m = mean(CoefErr_RMS);
s = std(CoefErr_RMS);
xAxis = linspace(1, 6, 6);
xAxis(1) = 0;
figure()
semilogy(xAxis, CoefErr_RMS, 'Marker','o', ...
    'LineStyle', 'none', 'Color','r', 'MarkerFaceColor','r')
hold on;
semilogy(xAxis, truth, 'k-', 'LineWidth', 1.5)
grid on;
title('RMS error compared with truth. MC simualtion = ' + string(MC))
xlabel('Harmonics degree [n]')
ylabel('coefficient error')

figure()
semilogy(xAxis, truth, 'k-', 'LineWidth', 1.5)
hold on;
errorbar(xAxis, m , s, 'LineWidth', 1.5, 'Color', 'r', 'Marker', 'square')
set(gca,'YScale','log')
grid on;
title('RMS error compared with truth. MC simualtion = ' + string(MC))
xlabel('Harmonics degree [n]')
ylabel('coefficient error')



% plot RMS error vs RMS sigma
figure()
plot(xAxis, CoefErr_RMS, 'Marker','o', ...
    'LineStyle', 'none', 'Color','r', 'MarkerFaceColor','r')
hold all;
plot(xAxis, boundRMS_up, 'Color','g', 'LineWidth', 1.5)
grid on;
title('RMS error + 1\sigma bounding. MC simualtion = ' + string(MC))
xlabel('Harmonics degree [n]')
ylabel('coefficient error')

% plot error vs s
figure()
subplot(2, 1, 1)
plot(linspace(1, 26, 26), CoefErr(:,1:26), 'Marker','o', ...
    'LineStyle', 'none', 'Color','r', 'MarkerFaceColor','r')
hold all;
plot(linspace(1, 26, 26), bound_up(:, 1:26), 'Color','g', 'LineWidth', 1.5)
plot(linspace(1, 26, 26), bound_down(:, 1:26), 'Color','g', 'LineWidth', 1.5)
plot(linspace(1, 26, 26), boundC_up(:, 1:26), 'Color','b', 'LineWidth', 1.5)
plot(linspace(1, 26, 26), boundC_down(:, 1:26), 'Color','b', 'LineWidth', 1.5)
grid on;
title('C_{nm} Error + 1\sigma bounding. MC simualtion = ' + string(MC))
xlabel('Harmonics degree [n]')
ylabel('coefficient error')
legend('MC values', '\sigma bound', 'sigma C')

subplot(2, 1, 2)
plot(linspace(1, 20, 20), CoefErr(:,27:end), 'Marker','o', ...
    'LineStyle', 'none', 'Color','r', 'MarkerFaceColor','r')
hold all;
plot(linspace(1, 20, 20), bound_up(:,27:end), 'Color','g', 'LineWidth', 1.5)
plot(linspace(1, 20, 20), bound_down(:,27:end), 'Color','g', 'LineWidth', 1.5)
plot(linspace(1, 20, 20), boundC_up(:,27:end), 'Color','b', 'LineWidth', 1.5)
plot(linspace(1, 20, 20), boundC_down(:,27:end), 'Color','b', 'LineWidth', 1.5)
grid on;
title('S_{nm} Error + 1\sigma bounding. MC simualtion = ' + string(MC))
xlabel('Harmonics degree [n]')
ylabel('coefficient error')



%%                      FUNCTIONS
function [CoefErr] = computeRMS_coefErr(n_max, Nc, Ns, err)
    % Description: compute the RMS of the coefficient error
    CoefErr = zeros(1, n_max);
    CoefErr(1) = sqrt((err(1))^2);
    m  = 0;
    n = 2;
    for j =2:Nc
        CoefErr(n) = CoefErr(n) + (err(j))^2;
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
        CoefErr(n) = CoefErr(n)  + (err(j))^2;
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