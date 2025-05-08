close all;
clear;
clc;


% Constants (arbitrary units)
GM = 10;          % Gravitational parameter of big sphere
Gm = 5;           % Gravitational parameter of internal mass concentration
d = 3;            % x-displacement of the internal small mass
n_points = 100;   % Number of points along the orbit
R_orbit = 6;      % Radius of the circular orbit
R_big = 5;        % Radius of big sphere
R_small = 0.5;    % Radius of internal mass

% Orbit positions
theta = linspace(0, 2*pi, n_points);
x = R_orbit * cos(theta);
y = R_orbit * sin(theta);

% Storage for eigenvectors
eig_vecs = zeros(2, n_points);

for i = 1:n_points
    % Position vector
    xi = x(i);
    yi = y(i);

    % Big mass contribution
    r = sqrt(xi^2 + yi^2);
    V_M = (GM / r^5) * [3*xi^2 - r^2, 3*xi*yi; 3*xi*yi, 3*yi^2 - r^2];

    % Small internal mass contribution
    x_d = xi - d;
    r_d = sqrt(x_d^2 + yi^2);
    V_m = (Gm / r_d^5) * [3*x_d^2 - r_d^2, 3*x_d*yi; 3*x_d*yi, 3*yi^2 - r_d^2];

    % total gravity potential  
    V_T = ((GM+Gm) / r^5) * [3*xi^2 - r^2, 3*xi*yi; 3*xi*yi, 3*yi^2 - r^2];

    % Total gravity gradient tensor
    V = V_M + V_m;

    % Eigen decomposition
    [V_eigvec, V_eigval] = eig(V);
    [~, idx] = max(abs(diag(V_eigval)));  % Principal direction
    eig_vecs(:, i) = V_eigvec(:, idx);
end

% Plotting
figure;
hold on;
axis equal;

% Plot orbit path
plot(x, y, 'k--', 'LineWidth', 1);

% Plot radial lines (optional: every 10th point for clarity)
for i = 1:10:n_points
    plot([0, x(i)], [0, y(i)], 'b:');
end

% Plot eigenvectors
quiver(x, y, eig_vecs(1, :), eig_vecs(2, :), 0.5, 'r');

% Plot big sphere boundary
theta_circ = linspace(0, 2*pi, 200);
plot(R_big * cos(theta_circ), R_big * sin(theta_circ), 'k-', 'LineWidth', 2);

% Plot small internal mass boundary
plot(d + R_small * cos(theta_circ), R_small * sin(theta_circ), 'r-', 'LineWidth', 2);

% Plot centers
plot(0, 0, 'ko', 'MarkerSize', 10, 'DisplayName', 'Big Sphere Center');
plot(d, 0, 'ro', 'MarkerSize', 8, 'DisplayName', 'Hidden Mass');

xlabel('x'); ylabel('y');
title('Gravity Gradient Eigenvectors with Radial Lines and Sphere Boundaries');
legend show;
grid on;

