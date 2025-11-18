clear;
clc;
close all;
format long g;
addpath('../functions/')
% % addpath('../../QGG_gravEstim/src/')
% % addpath('../../QGG_navigation/data/')
addpath('../data/')
set(0,'defaultAxesFontSize',16);

%%                          NSM THRESHOLD
% Description: Compute the NSM vs LS methods threshold for attitude or
% position.
% Author: Sergio Coll Ibars
% Date: 04/25/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = "attitude";     % option: position / attitude

% % Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % % % Eros parameters
% % path = "HARMCOEFS_EROS_CD_1.txt";
% % name = "EROS";
% % [Cnm, Snm, Re] = readCoeff(path);
% % n_max  = 10;
% % normalized = 1;
% % GM =  459604.431484721;          % Point mass value    [m^3/s^2]
% % W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
% % W0 = 0;                          % Initial asteroid longitude
% % RA = deg2rad(11.363);            % Right Ascension     [rad]
% % DEC = deg2rad(17.232);           % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

% Initial conditions
Nthrs = 10;
radius = linspace(6*Re, 1.4*Re, Nthrs);
meanThrs = ones(n_max-1, Nthrs) * NaN;
for k = 1:length(radius)
    r = radius(k);
    % % r       = 1E3;               % [m]
    r0      = [r;0;0];           % [ACI]
    v0      = [0;0;sqrt(GM/r)];  % [ACI]
    
    % time vector
    n = sqrt(GM / r^3);    % Mean motion         [rad/s]
    T = (2 * pi / n);
    rev = 3;
    f = 1/10;
    t = linspace(0, rev*T, rev*T * f);
    Nt = length(t);
    
    % measurement uncertianty
    sigma = 1E-12;                          % [1/s^2]
    noise0 = zeros(9, Nt);
    
    % Integrate trajectory
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    Nx = 6;
    PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
    [~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
        W0, W, RA, DEC, 0), t, [r0;v0;PHI0], options);
    rn = state_t(:, 1:3)';
    vn = state_t(:, 4:6)';
    
    % Consider covariance
    Ar = 1;                           % [-]
    Pxc = zeros(Ncs-1, 3); Pc = zeros(3, 3);
    Pxc_NSM = zeros(Ncs - 1, 6); Pc_NSM = zeros(6, 6);
    for j = 1:3
        Pc(j, j) = (Ar)^2;
    end
    for j = 1:length(Pc_NSM)
        Pc_NSM(j, j) = Ar^4;
    end
    c_NSM = ones(6, 1).*Ar^2;
    
    sigma_n = 1E3 * ones(1, n_max);
    [~, Pp] = perturb_coeff(sigma_n, n_max, X);
    P0 = Pp(2:end, 2:end); 
    
    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    [~, Mxc_NSM, Mcc_NSM] = get_considerCov_apriori(P0, Pc_NSM, Pxc_NSM);
    Ax = 0;  Ax_NSM = 0;
    R0 = diag([sigma, sigma, sigma, sigma, sigma].^2);
    for j = 1:Nt
        fprintf('Loading ... %.2f\n % ', j/Nt * 100);
        % current position
        rn_ACI = rn(:, j);
        vn     = rn.*0;
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % computed meas. partials
        [~, Hx_ACI, ~] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                noise0, Cnm, Snm);
        hx = [Hx_ACI(1, 2:end); Hx_ACI(2, 2:end); Hx_ACI(3, 2:end);Hx_ACI(5, 2:end);...
            Hx_ACI(6, 2:end)];
    
        % compute Consider Params partials for LS and NSM 
        [Hc] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, ACAF_ACI);
        hpos = [Hc(1, :);Hc(2,:);Hc(3,:);Hc(5, :);Hc(6, :)];
        hap = compute_posPartials_2ndOrder(GM, rn_ACI(1), rn_ACI(2), rn_ACI(3));
        
        if(threshold == "attitude")
            % compute attitude partials. % Inertially fixed
            [Hrot_grad] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI, ACAF_ACI);
            Hrot = Hrot_grad;
            hrot = [Hrot(1, :);Hrot(2,:);Hrot(3,:);Hrot(5, :);Hrot(6, :)];
            hc = hrot;
        else
             hc = hpos; % consider parameters matrix. Position NSM
        end
    
        % LS covariance
        Ax  = Ax  + (hx' * inv(R0) * hx);
        Mxc = Mxc + (hx' * inv(R0) * hc);
        Mcc = Mcc + (hc' * inv(R0) * hc); 
    
        % NSM covariance
        C = null(hc');
        hx_NSM = C' * hx;
        hap_NSM = C' * hap;
        r  = C' * R0 * C;
        
        Ax_NSM = Ax_NSM + hx_NSM' * inv(r) * hx_NSM;
        Mxc_NSM = Mxc_NSM + (hx_NSM' * inv(r) * hap_NSM);
        Mcc_NSM = Mcc_NSM + (hap_NSM' * inv(r) * hap_NSM); 
    end
    if(rank(Ax) ~= Ncs -1)
        disp('NOT ALL STATES OBSERVABLE');
    end
    % compute final covariance at epoch time. LS
    Px = inv(Ax);
    Sxc = -Px * Mxc;
    
    % compute final covariance at epoch time. NSM
    Px_NSM = inv(Ax_NSM);
    Sxc_NSM = -Px_NSM * Mxc_NSM;
    
    % compute analytical sigma pos value
    C = diag(Px_NSM - Px);
    A = diag(Sxc_NSM * Sxc_NSM');
    B = diag(Sxc * Sxc');
    sigmaTh1 = A.*0;
    for j = 1:length(A)
        c = C(j);
        b = B(j);
        a = A(j);
       
        % approximation including only 'b'
        if(threshold == "position")
% %             sigmaTh1(j) = (1/sigma)*sqrt(c/b) * 1E-9;                    % [m / E]
            sigmaTh1(j) = sqrt(c/b);                                         % [m]
        else
% %             sigmaTh1(j) = (1/sigma)*sqrt(c/b) * 1E-9 * 3600 * 180 / pi;  % [arcSec / E]
            sigmaTh1(j) = sqrt(c/b) * 3600 * 180 / pi;                       % [arcSec]
        end
    end
    
    % rms form
    RMS_thrs  = computeRMS_coeffErr(n_max, Nc, Ns, [0;sigmaTh1], Cnm.*0, Snm.*0); 
    
    % matrix form
    [C_thrs, S_thrs] = list2mat(n_max, Nc, Ns, [0; sigmaTh1]);
    C_thrs(1, :) = C_thrs(1, :).*NaN; C_thrs(2, :) = C_thrs(2, :).*NaN;
    S_thrs(:, 1) = S_thrs(:, 1).*NaN; S_thrs(1:2, :) = S_thrs(1:2, :).*NaN;

    % compute mean
    meanThrs(:, k) = RMS_thrs(2:end)';
end

p = polyfit(radius, log10(mean(meanThrs)), 4);
log_y_fit = polyval(p, radius);
y_fit = 10.^(log_y_fit);

% Plot
figure;
hold all;
x = radius; y = log10(y_fit);
% Fill below the curve (from ymin to y)
y_fill_bottom = [repmat(min(y)-2, size(y)) fliplr(y)];
fill([x fliplr(x)], y_fill_bottom, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% Fill above the curve (from y to ymax)
y_fill_top = [y fliplr(repmat(max(y)+2, size(y)))];
fill([x fliplr(x)], y_fill_top, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none')

plot(radius, log10(y_fit), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);

colors = lines(n_max -1); % 'lines' gives visually distinct colors
for i = 1:n_max-1
    plot(radius, log10(meanThrs(i,:)), '-s', ... % -s means line with square marker
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineStyle', 'none', ... % no connecting lines
        'DisplayName', ['degree = ' num2str(i+1)]);
end

% Formatting
yticks = -2:1:4; % Example: from 10^0 to 10^3
set(gca, 'YTick', yticks);
yticklabels = arrayfun(@(v) ['10^{' num2str(v) '}'], yticks, 'UniformOutput', false);
set(gca, 'YTickLabel', yticklabels);

xlabel('orbit radius [m]');
if(threshold == "position"), ylabel('[m]'); else, ylabel('[arcseconds]'); end
xlim([1.4*Re, 6*Re]);
legend('Location','best');
grid on;
hold off;
if(threshold == "position"), title('Position error threshold'); else,  ...
     title('Attitude error threshold'); end