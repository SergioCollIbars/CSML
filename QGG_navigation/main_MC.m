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
email = 0;  % send plots by email 1 = yes / 0 = no

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

Mc       = 10;
error_MC = ones(20, 1E4, Mc) * NaN; 
cov_MC   = ones(1E4, 20, Mc) * NaN;

% run MC
for mc = 1:Mc
    disp('Monte Carlo Iteration ... ' +string(mc))
    % run main code
    run("main_nav.m");
    
    % take error & covariance
    error_iter = state(:, 1:16)' - X;
    Ns = length(X(:, 1));
    Nt = length(time);

    error_MC(1:Ns, 1:Nt, mc) = error_iter;
    cov_MC(1:Nt, 1:Ns*Ns, mc)   = P; 

    clearvars -except mc Mc error_MC cov_MC planetParams Ns Nt TIME ...
        email diaryFile
end

% convert to UTC time
jd = 2451545 + TIME / planetParams(3) / 86400;
humanReadableTime = datetime(jd, 'ConvertFrom', ...
    'juliandate');
humanReadableTime.Format = 'MMM dd, yyyy';
date_init = string(humanReadableTime(1));
date_end  = string(humanReadableTime(end));
humanReadableTime.Format = 'MMM dd';

time = humanReadableTime';

%% PLOTER
lw = 2;
color1 = [204, 0, 204]./256;     % violet
color2 = 'r';                    % red

% plot position results
scale = planetParams(2)/ 1E3;               % [km]
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
title('Position norm')
grid on;
ylim([1E-5, 1]);

% plot velocity results
scale = planetParams(3) * planetParams(2);   % [m/s]
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
title('Velocity norm')
grid on;

% plot SRP scaling factor
scale = 1;
figure();
for mc = 1:Mc
    cov_iter = cov_MC(1:Nt, 1:Ns*Ns, mc);
    s = zeros(Ns, Nt);
    for k = 1:Nt
        p        = reshape(cov_iter(k, :), [Ns, Ns]); 
        s(:, k)  = sqrt(diag(p));
    end
    error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
    semilogy(time, abs(error(7, :)).* scale, 'LineWidth', lw, ...
        'Color', color1);
    hold on;
    if(mc == 1)
        semilogy(time, 3.*s(7, :).* scale, 'LineWidth', lw, ...
            'Color', 'k');
        hold on;
    end
end
ylabel('[-]')
title('SRP scaling factor norm')
grid on;


% plot gradiometer bias
ttlabel = ["xx", "xy", "xz", "yy", "yz", "zz"];
scale = planetParams(3)/1E-9;               % [Eotvos]
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
        semilogy(time, abs(error(7+k, :)).* scale, 'LineWidth', lw, ...
            'Color', color1);
        hold on;
        if(mc == 1)
            semilogy(time, 3.*s(7+k, :).* scale, 'LineWidth', lw, ...
                'Color', 'k');
            hold on;
            title('B_{' + ttlabel(k) + '}');
            ylabel('[Eotvos]')
            grid on;
        end
    end
end
sgtitle('Gradiometer Bias')

% plot gradiometer bias
ttlabel = ["x", "y", "z"];
scale = (planetParams(3)^2) * planetParams(2);               % [m/s^2]
figure();
for mc = 1:Mc
    cov_iter = cov_MC(1:Nt, 1:Ns*Ns, mc);
    s = zeros(Ns, Nt);
    for k = 1:Nt
        p        = reshape(cov_iter(k, :), [Ns, Ns]); 
        s(:, k)  = sqrt(diag(p));
    end
    error    = squeeze(error_MC(1:Ns, 1:Nt, mc));
    for k = 1:3
        subplot(1, 3, k)
        semilogy(time, abs(error(13+k, :)).* scale, 'LineWidth', lw, ...
            'Color', color1);
        hold on;
        if(mc == 1)
            semilogy(time, 3.*s(13+k, :).* scale, 'LineWidth', lw, ...
                'Color', 'k');
            hold on;
            title('B_{' + ttlabel(k) + '}');
            ylabel('[m/s^2]')
            grid on;
        end
    end
end
sgtitle('Accelerometer Bias')


% send resutls through email
if(email)
    % load console output 
    consoleText = fileread(diaryFile);
    
    sendOpenFiguresByEmail(...
        'sergiocollibars@gmail.com', ...
        'Auto Report - MATLAB Figures', consoleText);
end
delete(diaryFile);  % clean up immediately