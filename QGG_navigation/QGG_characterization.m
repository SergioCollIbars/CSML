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
%%                      QGG CHARACTERIZATION
% Description: set of scripts to characterize the sensor in a specific
% orbit. It testes things such as: 1/f freq threshold or batch
% initialization. 


%%                  BATCH INITIALIZATION
% WARNING: from QGG_navigation.m need to comment the following:
%          clear , clc, noiseSeed file, xnot, tmax and cut meas
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

[planetParams, poleParams, Cmat, Smat, ~] = load_universe("CR3BP", ...
    [0, 1], 1/60);

% clear kernels
cspice_kclear

scaleP = planetParams(2)/1000;              % [km]
scaleV = planetParams(3)*planetParams(2);   % [m/s]
scaleM = planetParams(3)^2;                 % [1/s^2]

time_array = (1:10).*86400; % [sec]
noiseIter = 100;
tmin = 0.295881248026939;                    % [rad]

sigmaPos = 1E5/planetParams(2);                   % [m]
sigmaVel = 10/(planetParams(2)*planetParams(3));  % [m/s]

% output values
RMS_iters = zeros(noiseIter * 6 ,length(time_array));
postfit_iters = zeros(noiseIter * 6 ,length(time_array));
for j = 1:length(time_array)
    for k = 1:noiseIter
        disp('time iteration  ' +string(j) + '/' + string(length(time_array)))
        disp('noise iteration  ' +string(k) + '/' + string(noiseIter))
        tmax = time_array(j).*planetParams(3) + tmin;

        posE = unifrnd(-sigmaPos, sigmaPos, [3, 1]);
        velE = unifrnd(-sigmaVel, sigmaVel, [3, 1]);
        Xnot = [posE;velE];

        run("QGG_navigation.m")

        if isempty(warnMsg)
            % compute error
            err = state(:, 1:6)' - X;
    
            % compute RMS
            maxPos = 6*k;
            minPos = maxPos - 5;
            RMS_iters(minPos:maxPos, j)  = [rms(err(1:3, :) * scaleP, 2); ...
                rms(err(4:6, :) * scaleV, 2)];
    
            % compute RMS postfit
            maxPos = 6*k;
            minPos = maxPos - 5;
            postfit = posIter((6*count)-5:6*count, :).*scaleM;
            postfit_iters(minPos:maxPos, j) = rms(postfit, 2, "omitnan");
        else
            % compute RMS
            maxPos = 6*k;
            minPos = maxPos - 5;
            RMS_iters(minPos:maxPos, j)  = ones(6, 1) * NaN;
    
            % compute RMS postfit
            maxPos = 6*k;
            minPos = maxPos - 5;
            postfit_iters(minPos:maxPos, j) = ones(6, 1) * NaN;
        end
        clear tmax;
    end
end

[mean_error, sigma_error] = computeStatistics(RMS_iters, 6);
[mean_postfit, sigma_postfit] = computeStatistics(postfit_iters, 6);

% plot RMS error. MC
line = sqrt(diag(R0))* scaleM;
tt = ["X position", "Y position", "Z position", "X velocity", "Y velocity", "Z velocity"];
plot_batchCharacteristics(time_array, RMS_iters, noiseIter, line*NaN, tt)
sgtitle('RMS state error changing noise seed. Solving LIS conditon')

% plot RMS error. Error bar
plot_errorBar(time_array, mean_error, sigma_error, tt);
sgtitle('RMS state error changing noise seed. Solving LIS conditon')

% plot RMS postfit. MC
tt = ["\Gamma xx", "\Gamma xy", "\Gamma xz", "\Gamma yy", "\Gamma yz", "\Gamma zz"];
plot_batchCharacteristics(time_array, postfit_iters, noiseIter, line, tt)
sgtitle('RMS postfit changing noise seed. Solving LIS conditon')

% plot RMS postfit. Error bar
plot_errorBar(time_array, mean_postfit, sigma_postfit, tt);
sgtitle('RMS state error changing noise seed. Solving LIS conditon')

%%                  LOST IN SPACE SOLVER 
% WARNING: from QGG_navigation.m need to comment the following:
%          clear , clc, noiseSeed file, Xnot
clear;
clc;
close all;
format long g;

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')
[planetParams, ~, ~, ~, TIME, ~] = load_universe("CR3BP", [0, 1/5*1.4968], 1/30);
cspice_kclear

Nt = round((TIME(end)-TIME(1))*(1/10/planetParams(3)) + 1);

% sigma meas
<<<<<<< HEAD
sigma = [1, 1/sqrt(2)] * 1E-12;                                             % [1/s^2]
=======
sigma = [1, 1/2] * 1E-12;                                             % [1/s^2]
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
sigma = sigma./(planetParams(3)^2);                                   % [-]
R_QGG = [sigma(1), sigma(2), sigma(2), sigma(1), sigma(2), sigma(1)]; % [-]

Ax = [1E8, 5E7, 1E7, 1E5, 1E3]./planetParams(2);                      % [m]
sigmaVel = 10/(planetParams(2)*planetParams(3));                      % [m/s]
Ax = 5E7./planetParams(2);                      % [m]
epsPos = 1E5/planetParams(2);                                         % [-]
epsVel = 1/(planetParams(2)*planetParams(3));                         % [-]
<<<<<<< HEAD
iters = 20;
=======
iters = 100;
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074

N0_ax    = ones(6, length(Ax)) * NaN;
H_warn   = ones(1, length(Ax)) * NaN;
sg       = ones(6, Nt, length(Ax)) * NaN;
mg       = ones(6, Nt, length(Ax)) * NaN;
posIterMC  = ones(3, Nt, iters) * NaN;
for k = 1:length(Ax)
    H = ones(6, iters) * NaN;
    s = ones(6, Nt, iters) * NaN;
    count_warn = 0;
    for j = 1:iters
        % show iteration number
        disp('Ax number: ' + string(k) + '/' + string(length(Ax)));
        disp('iter number: ' + string(j) + '/' + string(iters));

        % generate initial deviation error
        sigmaPos = Ax(k);
        posE = unifrnd(sigmaPos - epsPos, sigmaPos + epsPos, [3, 1]);
        velE = unifrnd(sigmaVel - epsVel, sigmaVel + epsVel, [3, 1]);
        Xnot = [posE;velE];
    
        % run code
        run("QGG_nav_main.m");

        % save estimated position
        posIterMC(:, :, j) = X(1:3, :);
        
        % normalize post fit data
        for h = 1:6
            posf(h, :) =  posf(h, :)./R_QGG(h); 
        end
        
        % save sigma values
        s(:, :, j) = posf;

        % compute chi square test 
        idx  = isnan(posf(1, :));
        if(nnz(idx) == length(TIME))
            H(:, j)    = 1;
        else
            % compute chi square test & check std and mean
            for h = 1:6
                idx  = isnan(posf(h, :));
                hval = chi2gof(posf(h, ~idx),'cdf',{@normcdf,0,1}, 'Alpha', 0.01);
                
                [warnMsg2, ~] = lastwarn;
                if ~isempty(warnMsg2) || ~isempty(warnMsg)
                   hval = 1;
                   count_warn = count_warn + 1;
                end
                warnMsg2 = ''; % clear errors related with H test

                H(h, j) = hval;
            end
        end
    end
    % count number of zeros in H == agree null hypothesis (gaussian)
    for h = 1:6
        N0_ax(h, k) = nnz(~H(h, :)); 
    end
    % save sg
    sg(:, :, k) = std(s, 0, 3);
    mg(:, :, k) = mean(s, 3);

    % save number of H test warnings
    H_warn(k) = count_warn;
end

% plot statistics
tt = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", ...
    "\Gamma_{yy}", "\Gamma_{yz}", "\Gamma_{zz}"];
plot_nullTest_stats(Ax.*planetParams(2)./1E3, N0_ax./iters, tt);
sgtitle('Null hypothesis percentage from \chi^2 test at different position errors');

% plot std smoothered for the different errors
figure()
mw = 50;
for k = 1:6
    subplot(2, 3, k)
    for j = 1:length(Ax)
        plot(t, movmean(sg(k, :, j), mw), 'LineWidth', 1)
        hold on;
    end
    xlabel('TIME [-]')
    title(tt(k))
end
sgtitle('Residuals standart deviation')

figure()
mw = 50;
for k = 1:6
    subplot(2, 3, k)
    for j = 1:length(Ax)
        plot(t, movmean(mg(k, :, j), mw), 'LineWidth', 1)
        hold on;
    end
    xlabel('TIME [-]')
    title(tt(k))
end
sgtitle('Residuals mean')
legend('\Delta_x = 1E8 m', '\Delta_x = 5E7', '\Delta_x = 1E7', ...
    '\Delta_x = 1E5', '\Delta_x = 1E3');


%%              FUNCTIONS

function [mean, sigma] = computeStatistics(val_iters, blckNum)
    Nt = length(val_iters(1, :));
    N  = length(val_iters(:, 1))/blckNum;
    
    % define output val
    mean = zeros(blckNum, Nt);
    count = zeros(blckNum, Nt);
  
    % compute mean
    
    for j = 1:N
        for k = 1:Nt
            maxPos = j * blckNum;
            minPos = maxPos - (blckNum - 1);
            val = val_iters(minPos:maxPos, k);
            
            if(~isnan(val))
                mean(:, k) = mean(:, k) + val;
                count(:, k) = count(:, k) + ones(6, 1);
            end
        end
    end
    mean = mean./count;

    % compute std
    count = zeros(blckNum, Nt);
    s = zeros(blckNum, Nt);
    for j = 1:N
        for k = 1:Nt
            maxPos = j * blckNum;
            minPos = maxPos - (blckNum - 1);
            val = val_iters(minPos:maxPos, k);
            if(~isnan(val))
                s(:, k) = s(:, k) + (val - mean(:, k)).^2;
                count(:, k) = count(:, k) + ones(6, 1);
            end
        end
    end
    count = count - ones(6, Nt);
    sigma  = sqrt(s./(count));
end

function [] = plot_batchCharacteristics(time_array, val_iters, itersNum, line, tt)
    figure()
    index = [1, 3, 5, 2, 4, 6];
    for g = 1:itersNum
        maxPos = 6*g;
        minPos = maxPos - 5;
        val = val_iters(minPos:maxPos, :);
        for k = 1:3
            subplot(3, 2, index(k))
            semilogy(time_array./86400, val(k, :), ...
                'Marker','o','MarkerFaceColor','r', ...
                'LineStyle','none');
            hold on;
            semilogy(time_array./86400, line(k)*ones(1, length(time_array)), '--', 'LineWidth', 2)
            grid on;
            if(g == 1)
                xlabel('TIME [days]')
                ylabel('[Km]')
                title(string(tt(k)))
            end
            hold on;
        end
        for k = 4:6
            subplot(3, 2, index(k))
            semilogy(time_array./86400, val(k, :),'Marker','o','MarkerFaceColor','r', ...
                'LineStyle','none');
            hold on;
            semilogy(time_array./86400, line(k)*ones(1, length(time_array)), '--', 'LineWidth', 2)
            grid on;
            if(g == 1)
                xlabel('TIME [days]')
                ylabel('[m/s]')
                title(string(tt(k)))
            end
            hold on;
        end 
    end
end

function [] = plot_errorBar(x_axis, y_axis, error, tt)
    figure()
    index = [1, 3, 5, 2, 4, 6];
    for k = 1:3
        subplot(3, 2, index(k))
        errorbar(x_axis./86400, y_axis(k, :), error(k, :), 'LineWidth', 2, ...
            'Color', 'r', 'Marker','o', 'MarkerFaceColor', 'b');
        grid on;
        xlabel('TIME [days]')
        ylabel('[Km]')
        set(gca,'YScale','log')
        title(tt(k))
    end
    for k = 4:6
        subplot(3, 2, index(k))
        errorbar(x_axis./86400, y_axis(k, :), error(k, :), 'LineWidth', 2, ...
            'Color', 'r', 'Marker','o', 'MarkerFaceColor', 'b');
        grid on;
        xlabel('TIME [days]')
        ylabel('[m/s]')
        set(gca,'YScale','log')
        title(tt(k))
    end 
end

function [] = plot_nullTest_stats(xAxis, yAxis, tt)
    figure()
    index = [1, 3, 5, 2, 4, 6];
    for h = 1:6
        subplot(3, 2, index(h))
        semilogx(xAxis, yAxis(h, :), 'sq --', ...
            'LineWidth', 2, 'Color', 'b');
        grid on;
        xlabel('position error [km]')
        ylabel("[-]")
        title(tt(h))
    end
end

function [] = plot_nullTest_stats_v2(xAxis, yAxis, tt)
    figure()
    semilogx(xAxis, yAxis, 'sq --', ...
        'LineWidth', 2, 'Color', 'b');
    grid on;
    xlabel('position error [km]')
    ylabel("[-]")
    title(tt)
end