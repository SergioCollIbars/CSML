clear;
clc;
close all;

%%        GRAVITY POTENTIAL PLOTTER FOR DIFFERENT ASTEROIDS

% PATHS
addpath('functions/');
addpath('data/');
addpath('../../QGG/data_files/')

% INPUTS
N = 500;
x = linspace(-500, 500, N);
y = linspace(-500, 500, N);
[X, Y] = meshgrid(x, y);

% PLANET
GM = 5.2;
Re = 246;
n_max = 6;
normalized = 0;
path = "HARMCOEFS_BENNU_OSIRIS_0.txt";
[C, S] = readCoeff(path);

% MATRICES
potential = ones(N, N) * NaN;
for j = 1:N % X direction
    for i = 1:N % Y direction
        r = [X(i, j); Y(i, j); 0];
        if(vecnorm(r) > Re)
            [~, dU, ~] = potentialGradient_nm(C, S, n_max, r, ...
                Re, GM, normalized);
            potential(i, j) = vecnorm(dU);
        end
    end
end


% PLOT
figure()
surf(X, Y, potential, 'EdgeColor','none');

