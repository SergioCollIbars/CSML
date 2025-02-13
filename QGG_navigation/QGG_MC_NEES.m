clear;
clc;
close all;
format long g;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
%%                    QGG NAVIGATION MC AND NEES TEST
% Description: compute the NEES test and plot the state error for a MC
% analysis.
% Author: Sergio Coll Ibars
% Date: 11/12/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nmc  = 25;
<<<<<<< HEAD
Nstate = 7;
Mc_error = ones(Nstate * Nmc, 1E7) * NaN;
=======
Mc_error = ones(6 * Nmc, 1E7) * NaN;
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
err_NEST = ones(Nmc, 1E7) * NaN;

% load initial covariance and generate inital errors
P0 = reshape(load('initCov_SRP.mat', 'p').p, [Nstate,Nstate]);
P0 = round(P0 * 1e24) / 1e24;
init_errors = mvnrnd(zeros(Nstate, 1), P0, Nmc);

for j = 1:Nmc
    disp('Iteration ' + string(j) + '/' + string(Nmc))
    
    % current initial state
    X0 = [load('initState_truth_SRP.mat').s;1.3] + init_errors(j, :)';

    % run code
    run("QGG_nav_main.m");
    
    % compute error
    err = state(:, 1:Nstate)' - X(1:Nstate, :);

    % save error
    maxVal = Nstate * j;
    minVal = maxVal - (Nstate - 1);
    Mc_error(minVal:maxVal, 1:length(TIME)) = err; 

    for kk = 1:length(TIME)
        p_k = reshape(P(kk, :), [Nstate,Nstate]);
        P_k = p_k(1:Nstate,1:Nstate);
        err_NEST(j, kk) = err(:, kk)' * inv(P_k) * err(:, kk);
    end
end

% compute sigma
cov  = zeros(Nstate, length(t));
cov2 = zeros(Nstate, length(t));
for j = 1:length(t)
    p = reshape(P(j, :), [Nstate,Nstate]);
    a = sqrt(diag(p));
    cov(:, j) = a(1:Nstate);
    cov2(:, j) = diag(p);
end

% plot results
lw1 = 1.5;
lw2 = 2;
<<<<<<< HEAD
color1 = [204, 0, 204]./256;     % violet
% % color1 = [0 0 1];
=======
% % color1 = [204, 0, 204]./256;     % violet
color1 = [0 0 1];
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
color2 = "#FF0000";              % red
color3 = "k";                    % black
set(0,'defaultAxesFontSize',16);

if(system == "EPHEM")
    jd = 2451545 + t / planetParams(3) / 86400;
    humanReadableTime = datetime(jd, 'ConvertFrom', ...
        'juliandate');
    humanReadableTime.Format = 'MMM dd, yyyy';
    date_init = string(humanReadableTime(1));
    date_end  = string(humanReadableTime(end));
    humanReadableTime.Format = 'MMM dd';

    time = humanReadableTime';
    xlb = "date";
    tt  = ' from ' + date_init + ' - ' + date_end;
else
    time = t'/planetParams(3)/86400;
    xlb = "days";
    tt = ' ';
end

% plot state error components
<<<<<<< HEAD
figure()
=======
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
index = [1, 3, 5];
scale = planetParams(2)./1000;
for j = 1:3
    for k = 1:Nmc
        maxVal = Nstate * k;
        minVal = maxVal - (Nstate - 1);
        err = Mc_error(minVal:maxVal, 1:length(TIME));

        subplot(3, 2, index(j))
        plot(time, err(j, :) * scale, 'LineWidth', lw1, 'Color', color1)
        merr = 0;
        hold all;
    end
    upper_bound = +3*cov(j, :) * scale;
    lower_bound = -3*cov(j, :) * scale;
    fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
            color1, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    xlabel('date')
    ylabel('error R_' + string(j) + '[km]')
    
    d = rms(err(j, :) * scale);
    d = round(d,2,"significant");
    str = {'RMS = ' + string(d)};
    legend(str,'Location', 'best');
end

index = [2, 4, 6];
scale = planetParams(3) * planetParams(2);
for j = 1:3
    for k = 1:Nmc
        maxVal = Nstate * k;
        minVal = maxVal - (Nstate - 1);
        err = Mc_error(minVal:maxVal, 1:length(TIME));

        subplot(3, 2, index(j))
        plot(time, err(j+3, :) * scale, 'LineWidth', lw1, 'Color', color1)
        merr = 0;
        hold all;
    end
    upper_bound = +3*cov(j+3, :) * scale;
    lower_bound = -3*cov(j+3, :) * scale;
    fill([time, fliplr(time)], [upper_bound, fliplr(lower_bound)], ...
            color1, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    xlabel('date')
    ylabel('error V_' + string(j) + '[m/s]')

    d = rms(err(j+3, :) * scale);
    d = round(d,2,"significant");
    str = {'RMS = ' + string(d)};
    legend(str,'Location', 'best');
end
sgtitle("Estimate state error + 3\sigma bounds over time")

% plot error magnitude
figure()
subplot(1, 2, 1)
d = sqrt(sum(cov2(1:3, :), 1));
scale = planetParams(2)./1000;
<<<<<<< HEAD
semilogy(time, 3.*d.* scale, 'LineWidth', lw2, 'Color', color3);
hold on;
for k = 1:Nmc
    maxVal = Nstate * k;
    minVal = maxVal - (Nstate - 1);
    err = Mc_error(minVal:maxVal, 1:length(TIME));

    semilogy(time, vecnorm(err(1:3, :)).* scale, 'LineWidth', lw1, 'Color',...
        color1)
    hold on;
end
=======
for k = 1:Nmc
    maxVal = 6 * k;
    minVal = maxVal - 5;
    err = Mc_error(minVal:maxVal, 1:length(TIME));

    semilogy(time, vecnorm(err(1:3, :)).* scale, 'LineWidth', lw, 'Color',...
        color1)
    hold all;
end
semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
grid on;
ylabel('[Km]')
title('Position norm')

subplot(1, 2, 2)
d = sqrt(sum(cov2(4:6, :), 1));
scale = planetParams(3) * planetParams(2);
<<<<<<< HEAD
semilogy(time, 3.*d.* scale, 'LineWidth', lw2, 'Color', color3);
hold on;
for k = 1:Nmc
    maxVal = Nstate * k;
    minVal = maxVal - (Nstate - 1);
    err = Mc_error(minVal:maxVal, 1:length(TIME));

    semilogy(time, vecnorm(err(4:6, :)).* scale, 'LineWidth', lw1, 'Color',...
        color1)
    hold on;
end
grid on;
ylabel('[m/s]')
title('Velocity norm')
legend('3 \sigma', 'error')
sgtitle('State error vector norm and 3 \sigma bound')

% plot augmented state
if(Nstate == 7)
    figure()
    semilogy(time, 3.*cov(7, :), 'LineWidth', lw2, 'Color', color3);
    hold on;
    for k = 1:Nmc
        maxVal = Nstate * k;
        minVal = maxVal - (Nstate - 1);
        err = Mc_error(minVal:maxVal, 1:length(TIME));

        semilogy(time, err(7, :), 'LineWidth', lw1, 'Color', color1)
        hold on;
    end
    grid on;
    xlabel('date')
    ylabel('[]')
    title('SRP \eta error and 3\sigma bound')
    legend('\eta error', '3\sigma')
end
=======
for k = 1:Nmc
    maxVal = 6 * k;
    minVal = maxVal - 5;
    err = Mc_error(minVal:maxVal, 1:length(TIME));

    semilogy(time, vecnorm(err(4:6, :)).* scale, 'LineWidth', lw, 'Color',...
        color1)
    hold all;
end
semilogy(time, 3.*d.* scale, 'LineWidth', lw, 'Color', 'k')
grid on;
ylabel('[m/s]')
title('Velocity norm')
legend('error', '3 \sigma')
sgtitle('State error vector norm and 3 \sigma bound')

>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
% save data
% % save('MC_state_err2.mat', 'Mc_error');
% % save('MC_state_cov2.mat', 'cov');

% NEST test
alpha = 0.05;                   % false alarm probability
n = length(err(:, 1));          % number of states

r1 = chi2inv(alpha/2, Nmc * n)./Nmc;    % upper limit
r2 = chi2inv(1-alpha/2, Nmc * n)./Nmc;  % lower limit

NEST  = mean(err_NEST(1:Nmc, 1:length(TIME)), 1);

figure()
plot(time, NEST(1:length(time)), 'o', 'MarkerSize', 2)
hold all;
plot(time, ones(1, length(time)).*r1, '--', 'LineWidth', 2)
plot(time, ones(1, length(time)).*r2, '--', 'LineWidth', 2)
title('NEES test \alpha = ' + string(alpha) + '. Mc = ' + string(Nmc));
xlabel('Time [days]')