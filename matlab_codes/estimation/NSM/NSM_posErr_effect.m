clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
addpath('../../../QGG_gravEstim/data_files/')
set(0,'defaultAxesFontSize',16);

%%            NSM POSITION ERROR EFFECT IN GRAV. ESTIMATION
% Description: Compute the trucation error for position errors depending on
% the order.
% Author: Sergio Coll
% Date: 10/09/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % path = "HARMCOEFS_EROS_CD_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % n_max  = 6;
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
r      = 0.3E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 0.5*[1;1;1];            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 1;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% noise values
noise0 = zeros(9, Nt);
sigma  = 1E-20;
noise  =  normrnd(0, sigma, [9, length(t)]);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;PHI0], options);
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;
vn = state_t(:, 4:6)';

% generate measurements
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm);

% Gravity estimation
Ax_R = 0; Ax_E = 0;
Nx_R = 0; Nx_E = 0;
R_R = diag([sigma, sigma, sigma, sigma, sigma].^2);
for j = 1:Nt
    % computed meas.
    [~, Hc_ACI, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);
    
    % RTN rotation matrix
    ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    rn_RTN = ACI_RTN' * rn(:, j);
    rn_ACI = rn(:, j);
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

     % compute Point mass approx error
    [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_ACI, ACAF_ACI);

    % Rummel's method
    [ax, nx] = rummels_method(Y(:, j), Hc_ACI, R_R, eye(3,3), noise(:, j));
    Ax_R  = Ax_R + ax;
    Nx_R  = Nx_R + nx;

    % error value
    hc = [Hc_ACI(1, 1:end); Hc_ACI(2, 1:end); Hc_ACI(3, 1:end);Hc_ACI(5, 1:end);...
        Hc_ACI(6, 1:end)];
    E = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)] * Ar;

    Ax_E  = Ax_E + ax;
    Nx_E  = Nx_E + hc' * inv(R_R) * E;
end
Xhat    = Ax_R\Nx_R;
Ac_err  = Ax_E\Nx_E;
sigma = sqrt(diag(inv(Ax_R)));

% Plot trucation error
[num_C, num_S, str_C, str_S] = SH_xlabel(n_max); 

figure()
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'k')
hold all;
semilogy(1:Nc-1, abs(Xhat(2:Nc)), 'Marker','square', 'LineStyle','--', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r')
semilogy(1:Nc-1, sigma(2:Nc), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b')
title('SH estimation. Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;
legend('truth', 'C_* error', '1 \sigma')

figure()
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'k')
hold all;
semilogy(1:Nc-1, abs(X(2:Nc) - Xhat(2:Nc)), 'Marker','square', 'LineStyle','--', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r')
semilogy(1:Nc-1, abs(Ac_err(2:Nc)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b')
title('Estimation error. Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;
legend('truth', 'Est. error', 'Pred. error')

%%  ESTIMATION ERROR @ DIFFERENT RADIUS
Nr     = 50;                             % radius points 
Nt     = 1E3;                            % points along the orbit
r      = linspace(Re*1.14, 5*Re, Nr);    % orbit radius vector
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
Ar = 0.5.*[1;1;1];
At = 1E-4.*[1;1;1];

RMSerr = ones(n_max, Nr) * NaN;
n2_err = ones(4, Nr) * NaN;
n3_err = ones(7,Nr) * NaN;
x     = ones(Nr, Nt) * NaN; 
y     = ones(Nr, Nt) * NaN; 
z     = ones(Nr, Nt) * NaN; 
for k = 1:Nr
    disp('Orbit ' +string(k) + '/' +string(Nr));

    % generate circular orbit
    rk = r(k);                                  % current radius [m]
    n = sqrt(GM/rk^3);                          % mean motion [rad/s]
    T = (2 * pi / n);                           % orbit period [s]
    t = linspace(0, rev*T, Nt);
    theta  = n.*t;                              % mean anomaly [rad]
    pos = [zeros(1, Nt); rk.*cos(theta); rk.*sin(theta)];  % pos. ACI coordinates

    % save orbit
    x(k, :) = pos(1, :);
    y(k, :) = pos(2, :);
    z(k, :) = pos(3, :);

% %     % compute position error effect in the SH estimation
% %     [Xerr] = compute_error_pos(t, pos, Nt, R_R, asterParams, poleParams, Cnm, Snm, Ar);

    % compute attitude error effect in the SH estimation
    [Xerr] = compute_error_att(t, pos, Nt, R_R, asterParams, poleParams, Cnm, Snm, At);

    % compute RMS error
    RMSerr(:, k) = computeRMS_coeffErr(n_max, Nc, Ns, ...
        -Xerr, Cnm.*0, Snm.*0); 
    n2_err(:, k) = abs([Xerr(2:4);Xerr(Nc+2)]);
    n3_err(:, k) = abs([Xerr(5:8);Xerr(Nc+3:Nc+5)]);
end

RMSt = computeRMS_coeffErr(n_max, Nc, Ns, X.*0, Cnm, Snm).*ones(Nr, n_max);
n2_true = abs([X(2:4);X(Nc+2)]).*ones(4, Nr);
n3_true = abs([X(5:8);X(Nc+3:Nc+5)]).*ones(7, Nr);

n2err_perct = abs(n2_err./n2_true.*100);
n3err_perct = abs(n3_err./n3_true.*100);

%% PLOT. LINES

% moving mean window
c = 'r';
lw = 2;

% RMS error
figure()
for n = 1:n_max
    subplot(2, 3, n)
    semilogy(r./1E3, smoothdata(RMSerr(n, :)), 'LineWidth', lw, 'Color', c)
    hold on; 
    semilogy(r./1E3, RMSt(:, n)', 'LineWidth', lw, 'Color', 'k')
    xlabel('radius [km]')
    ylabel('RMS error')
    title('n = ' + string(n))
end
sgtitle('Error introduced by ' + string(vecnorm(Ar)) + ' m position error for each SH RMS degree')
legend('error', 'truth')

% J2 error
tt2 = ["C_{20}", "C_{21}", "C_{22}", "S_{21}"];
figure()
for n = 1:4
    subplot(2, 2, n)
    plot(r./1E3,  smoothdata(n2_err(n, :)), 'LineWidth', lw, 'Color', c)
    hold on; 
    plot(r./1E3, n2_true(n, :), 'LineWidth', lw, 'Color', 'k')
    xlabel('radius [km]')
    ylabel('Error')
    title(tt2(n))
end
sgtitle('Error introduced by ' + string(vecnorm(Ar)) + ' m position error for degree 2')

% J3 error
tt3 = ["C_{30}", "C_{31}", "C_{32}", "C_{33}", "S_{31}", "S_{32}", "S_{33}"];
figure()
for n = 1:7
    subplot(2, 4, n)
    plot(r./1E3,  smoothdata(n3_err(n, :)), 'LineWidth', lw, 'Color', c)
    hold all; 
    plot(r./1E3, n3_err(n, :), 'LineWidth', 1, 'Color', 'b')
    plot(r./1E3, n3_true(n, :), 'LineWidth', lw, 'Color', 'k')
    xlabel('radius [km]')
    ylabel('Error')
    title(tt3(n))
end
sgtitle('Error introduced by ' + string(vecnorm(Ar)) + ' m position error degree 3')

%% FUNCTIONS
function [num_C, num_S, str_C, str_S] = SH_xlabel(n_max)    
    num_C = ones(1, n_max-1) * NaN;
    num_S = num_C;

    str_C = cell(1, n_max - 1);
    str_S = str_C;

    num_C(1) = 1;
    for j = 3:n_max
        num_C(j-1) = j + num_C(j-2);
    end
    
    num_S(1) = 1;
    for j = 3:n_max
        num_S(j-1) = (j-1) + num_S(j-2);
    end

    for j = 2:n_max
        str_C{j - 1} = "C_{" + string(j) + "0}";
        str_S{j - 1} = "S_{" + string(j) + string(j) + "}";
    end 
end


function [ax, nx] = rummels_method(Y, Hc, R, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;
    
    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];
    hc = [Hc(1, 1:end); Hc(2, 1:end); Hc(3, 1:end);Hc(5, 1:end);...
        Hc(6, 1:end)];
    
    % information and normal matrices
    ax = hc' * inv(R) * hc;
    nx = hc' * inv(R) * dY;
end


function [Xerr] = compute_error_pos(t, pos, Nt, R, asterParams, poleParams, Cnm, Snm, Ar)
    Ax = 0; Nx = 0;
    noise0 = zeros(9, 1);
    
    GM         = asterParams(1);
    Re         = asterParams(2);
    n_max      = asterParams(3);
    normalized = asterParams(4);
    
    W   = poleParams(1);
    W0  = poleParams(2);
    RA  = poleParams(3);
    DEC = poleParams(4);

    for j = 1:Nt
        % computed meas.
        [~, Hc_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [pos(:, j)', 0,0,0], ...
                noise0, Cnm, Snm);
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
         % compute Point mass approx error
        [Hp0] = compute_posPartials(0, normalized, Cnm, Snm, Re, GM, pos(:, j), ACAF_ACI);
        [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, pos(:, j), ACAF_ACI);
        % % Hp  = Hp - Hp0;

        % error value
        hc = [Hc_ACI(1, 1:end); Hc_ACI(2, 1:end); Hc_ACI(3, 1:end);Hc_ACI(5, 1:end);...
            Hc_ACI(6, 1:end)];
        E = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)] * Ar;
    
        Ax  = Ax + hc' * inv(R) * hc;
        Nx  = Nx + hc' * inv(R) * E;
    end
    Xerr = Ax\Nx;
end


function [Xerr] = compute_error_att(t, pos, Nt, R, asterParams, poleParams, Cnm, Snm, At)
    Ax = 0; Nx = 0;
    noise0 = zeros(9, 1);
    
    GM         = asterParams(1);
    Re         = asterParams(2);
    n_max      = asterParams(3);
    normalized = asterParams(4);
    
    W   = poleParams(1);
    W0  = poleParams(2);
    RA  = poleParams(3);
    DEC = poleParams(4);

    for j = 1:Nt
        % computed meas.
        [~, Hc_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [pos(:, j)', 0,0,0], ...
                noise0, Cnm, Snm);
        
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(j);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
         % compute Point mass approx error
        [Hrot] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, pos(:, j), ACAF_ACI);
    
        % error value
        hc = [Hc_ACI(1, 1:end); Hc_ACI(2, 1:end); Hc_ACI(3, 1:end);Hc_ACI(5, 1:end);...
            Hc_ACI(6, 1:end)];
        E = [Hrot(1, :);Hrot(2,:);Hrot(3,:);Hrot(5, :);Hrot(6, :)] * At;
    
        Ax  = Ax + hc' * inv(R) * hc;
        Nx  = Nx + hc' * inv(R) * E;
    end
    Xerr = Ax\Nx;
end


