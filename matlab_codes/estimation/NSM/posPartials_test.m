%%                     POSITION PARTIALS PARTIALS TEST
%   Description: Compute the error in gradiometer meas. for first and
%   second order position partials. 
%   Author: Sergio Coll Ibars
%   Date: 02/10/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% Initial conditions & position error
r      = 0.5E3;             % [m]
phi    = pi/2;              % [rad]
lambda = 0;                 % [rad]
theta  = pi/2 - phi;        % Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

Ar = 5*[1;1;1].*1;            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% Integrate trajectory
options = odeset('RelTol',1e-11,'AbsTol',1e-11);
STM0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;                                        % constant position error
% % rn = state_t(:, 1:3)' + [sin(1E-4.*t);sin(1E-3.*t);sin(5E-4.*t)].*Ar;       % sinusoidal position error
vn = state_t(:, 4:6)';

% generate measurements @ nominal and true trajectory. ACI frame
noise0 = zeros(9, Nt);
[Ytrue, ~, ~]    = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                    noise0, Cnm, Snm);
[Ynominal, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, [rn' ,vn'], ...
                    noise0, Cnm, Snm);

% compute the position error expansion
err1 = ones(5, Nt) * NaN;
err2 = ones(5, Nt) * NaN;
for j = 1:Nt
     % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % gradiometer error
    dY = [Ytrue(1:3, j); Ytrue(5:6, j)] - [Ynominal(1:3, j); Ynominal(5:6, j)];

    % compute position partials. 1st and 2nd order
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, j), ACAF_ACI);
    hp2 = compute_posPartials_2ndOrder(GM, rn(1, j), rn(2, j), rn(3, j));
    hp1 = [Hpos(1:3, :);Hpos(5:6, :)];
    
    % compute error
    T = Ar * Ar';
    d = [T(1,1);T(1,2);T(1,3);T(2,2);T(2,3)];
    err1(:, j)  = dY + (hp1 * Ar);
    err2(:, j)  = dY + (hp1 * Ar + hp2 * ones(6, 1).*d);
end

% plot error over time
colors = [1 0 0;     % Red
          0 1 0;     % Green
          0 0 1;     % Blue
          0.5 0.5 0; % Olive (custom RGB)
          1 0 1];    % Magenta
figure()
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren'); % Set color order
plot(t./3600, err1./1E-9, 'LineWidth', 2)
hold on;
plot(t./3600, err2./1E-9, 'LineWidth', 2, 'LineStyle', '--')
xlabel('Time [hours]');
ylabel('[Eotvos]')
legend('\Gamma_{xx}', '\Gamma_{xy}', '\Gamma_{xz}', '\Gamma_{yy}', '\Gamma_{yz}')
grid on;
title('Taylor expansion error for gradiometer measurements \Delta_r = ' + string(Ar(1)) + ' m')
xlim([10, 15]);

% extra text
xPos = 0.4; % X position in normalized figure coordinates
yPos = 0.8; % Y position in normalized figure coordinates
width = 0.1; % Width of the text box
height = 0.1; % Height of the text box

% Add annotation text box
annotation('textbox', [xPos, yPos, width, height], 'String', {'Solid line: O(2)', 'Dashed line: O(3)'}, ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
