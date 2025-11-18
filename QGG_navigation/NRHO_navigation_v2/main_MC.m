clear;
clc;
close all;
format long g;

% email settings
setpref('Internet','E_mail','sergiocollibars@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','sergiocollibars@gmail.com');
setpref('Internet','SMTP_Password','ptcq ybug jcgh wajg');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true'); % For TLS
props.setProperty('mail.smtp.port','587'); % Common TLS port
email = 1;  % send plots by email 1 = yes / 0 = no

% Start diary logging
diaryFile = fullfile(tempdir, 'console_log.txt');
if isfile(diaryFile)
    delete(diaryFile);
end
diary(diaryFile);
diary on;

% plot settings
set(0,'defaultAxesFontSize',16);

%%              RUN MAIN_NAV.M CODE IN A MC BATCH
% Description: generate MC samples for the main_nav.m file which runs
% navigation in cislunar space using Ephemerides files
%
% WARNING: Make sure plotResults = 0 and loadData = 1 before RUN
%
% Date: 10/01/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mc       = 1;
error_MC = ones(20, 1E4, Mc) * NaN; 
cov_MC   = ones(1E4, 20, Mc) * NaN;
Ns       = 12;

% run MC
for mc = 1:Mc
    disp('Monte Carlo Iteration ... ' +string(mc))
    % run main code
    run("main_nav.m");
    
    % take error & covariance
    error_iter = state(:, 1:Ns)' - X;
    Ns = length(X(:, 1));
    Nt = length(time);

    error_MC(1:Ns, 1:Nt, mc) = error_iter;
    cov_MC(1:Nt, 1:Ns*Ns, mc)   = P; 

    gamma_old = gamma;

    clearvars -except mc Mc error_MC cov_MC planetParams Ns Nt TIME ...
        email diaryFile r_sc_moon gamma_old
end

% compute NEES test
NEES = zeros(1, length(TIME));
NEES_pos_vel = zeros(1, length(TIME));
NEES_bias    = zeros(1, length(TIME));
alpha = 0.05;
for k = 1:length(TIME)
    ek = squeeze(error_MC(1:Ns, k, :));
    Pk = reshape(squeeze(cov_MC(k, :, 1)), [Ns, Ns]);

    ek_pos_vel = ek(1:6, :);
    Pk_pos_vel = Pk(1:6, 1:6);

    ek_bias    = ek(7:end, :);
    Pk_bias    = Pk(7:12, 7:12);
    for mc = 1:Mc
        NEES(k) = NEES(k) + ek(:, mc)' * (Pk \ ek(:, mc));
        NEES_pos_vel(k) = NEES_pos_vel(k) + ...
            ek_pos_vel(:, mc)' * (Pk_pos_vel \ ek_pos_vel(:, mc));
        NEES_bias(k) = NEES_bias(k) + ...
            ek_bias(:, mc)' * (Pk_bias \ ek_bias(:, mc));
    end
end
NEES = NEES./Mc;
NEES_pos_vel = NEES_pos_vel./Mc;
NEES_bias    = NEES_bias./Mc;
r1 = chi2inv(alpha/2, Mc * Ns)/Mc;
r2 = chi2inv(1-alpha/2, Mc * Ns)/Mc;

r1_pos_vel = chi2inv(alpha/2, Mc * 6)/Mc;
r2_pos_vel = chi2inv(1-alpha/2, Mc * 6)/Mc;

r1_bias = chi2inv(alpha/2, Mc * 6)/Mc;
r2_bias = chi2inv(1-alpha/2, Mc * 6)/Mc;

% % [gamma_new] = compute_gamma(NEES, Ns);
% % gamma = gamma_old.*gamma_new;
% % % save gamma value
% % save("gamma.mat", "gamma");

% convert to UTC time
jd = 2451545 + TIME / planetParams(3) / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', ...
    'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';

time = humanReadableTime';

% get apo-apsis and peri-apsis passages
idx_min = islocalmin(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations
idx_max = islocalmax(r_sc_moon, 'MinProminence', 0.0005);  % ignores small fluctuations


% save covariance file
F = squeeze(cov_MC(:, :, 1));
data = {F, time};
save('cov_Mc_1.mat', 'data');

%% PLOTER
lw = 2;
color1 = [204, 0, 204]./256;     % violet
color2 = 'r';                    % red

% Plot NEES test
figure();
plot(time, NEES, 'LineWidth', 2, 'Color','b')
hold all;
% % plot(time, gamma, 'LineWidth', 1.5, 'Color','k')
yline(r1, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
yline(r2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
grid on;
title('NEES test');

figure()
subplot(1, 2, 1)
plot(time, NEES_pos_vel, 'LineWidth', 2, 'Color','r')
hold all;
yline(r1_pos_vel, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
yline(r2_pos_vel, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
grid on;
title('NEES test pos & vel');

subplot(1, 2, 2)
plot(time, NEES_bias, 'LineWidth', 2, 'Color','r')
hold all;
yline(r1_bias, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
yline(r2_bias, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--')
grid on;
title('NEES test bias');


% plot position results
scale = planetParams(2) / 1E3;               % [km]
figure();
for mc = 1:Mc
    cov_iter = cov_MC(1:Nt, 1:Ns*Ns, mc);
    s = zeros(Ns, Nt);
    for k = 1:Nt
        p        = reshape(cov_iter(k, :), [Ns, Ns]); 
        s(:, k)  = sqrt(diag(p));
    end
    error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
    semilogy(time, vecnorm(error(1:3, :)).* scale, 'LineWidth', lw, ...
        'Color', color1);
    hold on;
    if(mc == 1)
        semilogy(time, 3.*vecnorm(s(1:3, :)).* scale, 'LineWidth', lw, ...
            'Color', 'k');
        hold on;
    end
end
ylabel('[Km]')
title('Position error norm & 3\sigma bound')
grid on;
ylim([1E-5, 1]);
if(sum(idx_min)~=0)
    xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
        'LineStyle', '--')
end
if(sum(idx_max)~=0)
    xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
        'LineStyle', '--')
end
% plot velocity results
scale = planetParams(2) * planetParams(3);   % [m/s]
figure();
for mc = 1:Mc
    cov_iter = cov_MC(1:Nt, 1:Ns*Ns, mc);
    s = zeros(Ns, Nt);
    for k = 1:Nt
        p        = reshape(cov_iter(k, :), [Ns, Ns]); 
        s(:, k)  = sqrt(diag(p));
    end
    error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
    semilogy(time, vecnorm(error(4:6, :)).* scale, 'LineWidth', lw, ...
        'Color', color1);
    hold on;
    if(mc == 1)
        semilogy(time, 3.*vecnorm(s(4:6, :)).* scale, 'LineWidth', lw, ...
            'Color', 'k');
        hold on;
    end
end
ylabel('[m/s]')
title('Velocity error norm & 3\sigma bound')
grid on;
if(sum(idx_min)~=0)
    xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
        'LineStyle', '--')
end
if(sum(idx_max)~=0)
    xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
        'LineStyle', '--')
end

% plot gradiometer bias
ttlabel = ["xx", "xy", "xz", "yy", "yz", "zz"];
scale = 1E6;               % [milli-Eotvos]
figure();
for mc = 1:Mc
    cov_iter = cov_MC(1:Nt, 1:Ns*Ns, mc);
    s = zeros(Ns, Nt);
    for k = 1:Nt
        p        = reshape(cov_iter(k, :), [Ns, Ns]); 
        s(:, k)  = sqrt(diag(p));
    end
    error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
    for k = 1:6
        subplot(3, 2, k)
        plot(time, error(6+k, :).* scale, 'LineWidth', 1.2, ...
            'Color', color2);
        hold on;
        if(mc == 1)
            plot(time, 3.*s(6+k, :).* scale, 'LineWidth', lw, ...
                'Color', 'k');
            hold on;
            plot(time, -3.*s(6+k, :).* scale, 'LineWidth', lw, ...
                'Color', 'k');
            title('B_{' + ttlabel(k) + '}');
            ylabel('[mE]')
            grid on;
        end
    end
end
sgtitle('Gradiometer bias error & 3\sigma bound')

%% plot pos & vel in cartesian components
% plot position results
scale = planetParams(2);               % [m]
figure();
tt = ["X", "Y", "Z"];
cov_iter = cov_MC(1:Nt, 1:Ns*Ns, Mc);
for k = 1:Nt
    p        = reshape(cov_iter(k, :), [Ns, Ns]); 
    s(:, k)  = sqrt(diag(p));
end
for j = 1:3
    subplot(1, 3, j)
    for mc = 1:Mc
        error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
        plot(time, error(j, :).* scale, 'LineWidth', lw, ...
            'Color', color1);
        hold all;
        if(mc == 1)
            upper_bound = +3*s(j, :) * scale;
            lower_bound = -3*s(j, :) * scale;
            fill([time', fliplr(time')], [upper_bound, fliplr(lower_bound)], ...
                color1, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            hold all;
            plot(time, 3.*(s(j, :)).* scale, 'LineWidth', 0.4, ...
            'Color', 'k');
            plot(time, -3.*(s(j, :)).* scale, 'LineWidth', 0.4, ...
            'Color', 'k');
        end
    end
    title(tt(j));
    ylabel('[m]')
    grid on;
    if(sum(idx_min)~=0)
        xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
            'LineStyle', '--')
    end
    if(sum(idx_max)~=0)
        xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
            'LineStyle', '--')
    end
end
sgtitle('Position error & 3\sigma bound')

% plot velocity results
scale = planetParams(2) * planetParams(3) * 1E3;   % [mm/s]
figure();
tt = ["X_x", "V_y", "V_z"];
for j = 1:3
    subplot(1, 3, j)
    for mc = 1:Mc
        error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
        plot(time, error(j+3, :).* scale, 'LineWidth', lw, ...
            'Color', color1);
        hold all;
        if(mc == 1)
            upper_bound = +3*s(j+3, :) * scale;
            lower_bound = -3*s(j+3, :) * scale;
            fill([time', fliplr(time')], [upper_bound, fliplr(lower_bound)], ...
                color1, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            hold all;
            plot(time, 3.*(s(j+3, :)).* scale, 'LineWidth', 0.4, ...
            'Color', 'k');
            plot(time, -3.*(s(j+3, :)).* scale, 'LineWidth', 0.4, ...
            'Color', 'k');
        end
    end
    title(tt(j));
    ylabel('[mm/s]')
    grid on;
    if(sum(idx_min)~=0)
        xline(time(idx_min), 'LineWidth', 1.1, 'Color', 'r', ...
            'LineStyle', '--')
    end
    if(sum(idx_max)~=0)
        xline(time(idx_max), 'LineWidth', 1.1, 'Color', 'b', ...
            'LineStyle', '--')
    end
end
sgtitle('Velocity error & 3\sigma bound')

%% send resutls through email
if(email)
    % load console output 
    consoleText = fileread(diaryFile);
    
    sendOpenFiguresByEmail(...
        'sergiocollibars@gmail.com', ...
        'Auto Report - MATLAB Figures', consoleText);
end
delete(diaryFile);  % clean up immediately

%% FUNCTION
function [gamma] = compute_gamma(NEES_k, n)
    % cosntants
    lambda = 0.98;
    d = 0.05;
    eta_up = 0.15;
    eta_down = 0.1;
    c_min = 0.6;
    c_max = 10;

    r_k = NEES_k./n;
    r_hat = zeros(1, length(NEES_k));
    gamma = ones(1, length(NEES_k));
    for k = 2:length(NEES_k)
        r_hat(k) = lambda*r_hat(k-1) + (1-lambda)*r_k(k);
        delta_k = r_hat(k) - 1;

        if(abs(delta_k) <=d)
            gamma(k) = 1; 
        elseif(delta_k > d)
            gamma(k) = exp(eta_up * delta_k); 
        elseif(delta_k < -d)
            gamma(k) = exp(eta_down * delta_k); 
        end

        gamma(k) = min([max([gamma(k), c_min]), c_max]);
    end
end