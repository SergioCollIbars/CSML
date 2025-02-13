clear;
close all;
clc;

%%          PARTIALS REPRESENTATION


% planet parameters
GM = 5.2;   % [m^3/s^2]
Re = 250;   % [m]

% create meshgrid
n = 100; % Define the number of points between a and b (adjust as needed)

x = linspace(Re, 6*Re, n);
y = linspace(GM, 6*GM, n);
[X, Y] = meshgrid(x, y);
Z = X.*0;

for j = 1:n
    for i = 1:n
        % current position value
        pos = X(i, j);
        gm  = Y(i, j);

        Z(i, j) = 2 * gm/pos^3;
    end
end

% plot figure
figure()
surf(X, Y, Z)
