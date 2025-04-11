clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              ODEST FIRST ERROR TRUNCATION ERROR
% Description: Compute the trucation error for position errors depending on
% the order.
% Author: Sergio Coll
% Date: 09/26/24

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

% Orbit coordinates.
N = 2;
r  = 0.6E3;       % Orbit radius     [m]
[P, L] = meshgrid(linspace(-pi/2, pi/2, N), linspace(0, 2*pi, N));

% position error
Ar = [1;1;1];   % [m]

err1 = ones(N, N) * NaN;
err2 = ones(N, N) * NaN;
for x = 1:N
    for y = 1:N
        phi   = P(x, y);
        lambda = L(x, y);
        theta  = pi/2 - phi;% Orbit colatitude [m]
        R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
            sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
            cos(theta), -sin(theta), 0];
        r0 = R * [r;0;0];
        v0 = R * [0;0;sqrt(GM/r)];
        
        % RTN rotation matrix
        ACI_RTN = RTN2ECI(r0, v0);
        rn_RTN = ACI_RTN' * r0; 
        
        % compute position partial. 1st order
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_RTN, eye(3,3));
        
        % compute Point mass approx error
        [Hpos0] = compute_posPartials(0, normalized, Cnm, Snm, Re, GM, rn_RTN, eye(3,3));
        
        E = [Hpos(5, :) - Hpos0(5, :); Hpos(4, :) - Hpos(6, :) -...
            Hpos0(4, :) + Hpos0(6, :)] * Ar;
        err1(x, y) = E(1);
        err2(x, y) = E(2);
    end
end

% Plot error
figure
surf(P, L, err1)
xlabel('Latitude, \phi')
ylabel('Longitude, \lambda')
zlabel('[1/s^2]')
title('Position sensitivity truncation error, radius = ' + string(r/1E3) + 'Km')
shading interp
view(0, 90)
colorbar

%%
% compute SH error associated with trucation
r      = 0.3E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]
Ar = 10*[1;1;1];            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 1;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0], options);
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;
vn = state_t(:, 4:6)';
plot_trajectory(state_t, 'BENNU');
figure()
plot(t./T, vecnorm(state_t(:, 1:3)'), 'LineWidth', 2)
hold on;
plot(t./T, ones(1, Nt)*Re, 'LineWidth', 2, 'Color', 'r', 'LineStyle','--')
xlabel('Orb. Period number, T')
ylabel('[m]')
title('Orbit radius along trajectory. T = ' + string(T./3600) + ' h')
legend('orbit radius', 'brillouin sphere')

Ax = 0;
Nx = 0;
noise0 = zeros(9, Nt);
sigma  = 6.32E-20;
R = diag([sigma, sigma].^2);
for j = 1:Nt
      % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % computed meas.
    [~, ~, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);
    hc = [Hc_RTN(8, 2:end); Hc_RTN(5, 2:end) - Hc_RTN(9, 2:end)];
    
    % RTN rotation matrix
    ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    rn_RTN = ACI_RTN' * rn(:, j);
    
    % compute position partial. 1st order
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);
    
    % compute Point mass approx error
    [Hpos0] = compute_posPartials(0, normalized, Cnm, Snm, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);
    
    E = [Hpos(8, :) - Hpos0(8, :); Hpos(5, :) - Hpos(9, :) -...
        Hpos0(5, :) + Hpos0(9, :)] * Ar;

    Ax  = Ax + hc' * inv(R) * hc;
    Nx  = Nx + hc' * inv(R) * E;
end
Ac_err = inv(Ax) * Nx;
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma = sqrt(diag(inv(Ax)));

% Plot trucation error
[num_C, num_S, str_C, str_S] = SH_xlabel(n_max); 
width  =  1130;
height = 503;

figure()
subplot(1, 2, 1)
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'k')
hold all;
semilogy(1:Nc-1, abs(Ac_err(1:Nc-1)), 'Marker','square', 'LineStyle','--', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r')
semilogy(1:Nc-1, sigma(1:Nc-1), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b')
title('Trucation error in SH estimation. Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;
legend('truth', 'C_* error', '1 \sigma')

% percentage error
subplot(1, 2, 2)
Ac_err(2) = 0;
bar(1:Nc-1, abs(Ac_err(1:Nc-1))./abs(X(2:Nc)).*100)
set(gca, 'YScale', 'log');
title('% error in SH estimation. Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;

pos = get(gcf, 'Position');  % Get current position: [left, bottom, width, height]
set(gcf, 'Position', [pos(1), pos(2), width, height]);  % Modify only width and height

figure()
subplot(1, 2, 1)
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'k')
hold on;
semilogy(1:Ns, abs(Ac_err(Nc:Nc+Ns-1)), 'Marker','square', 'LineStyle','--', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'r')
semilogy(1:Ns, sigma(Nc:Nc+Ns-1), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','b')
title('Trucation error in SH estimation. Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;

subplot(1, 2, 2)
bar(1:Ns, abs(Ac_err(Nc:Nc+Ns-1))./abs(X(Nc+1:Nc+Ns)).*100)
set(gca, 'YScale', 'log');
title('% error in SH estimation. Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;

pos = get(gcf, 'Position');  % Get current position: [left, bottom, width, height]
set(gcf, 'Position', [pos(1), pos(2), width, height]);  % Modify only width and height

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


