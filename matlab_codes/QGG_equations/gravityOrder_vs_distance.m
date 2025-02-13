clc;
close all;
clear;

%%                    GRAVITY ORDER VS DISTNACE
%
%   Description: This code computes the gravity influece of each
%   coefficient as a ratio vs the radial distance.


% IMPORTS
addpath('../../QGG_gravEstim/data_files/')
addpath('../../QGG_gravEstim/modules/estimation_module/functions')

% longitude meshgrid
lat = deg2rad(0);
lon = deg2rad(45);
N = 300;

% Planet parameters
Re = 246;                           % Planet radius [m]   6355E3 (polar radius)
GM = 5.2;                           % gravity param [m^3 / s^2]
H = linspace(Re, 20 * Re, N);     % SC altitude [m]

n_max = 6;                          % spherical harmonic order
normalized = 1;

path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[C_mat, S_mat, Re] = readCoeff(path);

% Gravity value
Un = zeros(6, length(H));

for k = 1:length(H)
    % get position vector
    x = H(k) * cos(lat) * cos(lon);
    y = H(k) * cos(lat) * sin(lon);
    z = H(k) * sin(lat);

    r = [x ; y; z];

    [~, dUf, ~] = potentialGradient_nm(C_mat, S_mat, 0, ...
                                                r, Re, GM, normalized);
    Un(1, k) = vecnorm(dUf);
end

for n = 2:6
    for k = 1:length(H)
        % get position vector
        x = H(k) * cos(lat) * cos(lon);
        y = H(k) * cos(lat) * sin(lon);
        z = H(k) * sin(lat);

        r = [x ; y; z];

        [~, dUf, ~] = potentialGradient_nm(C_mat, S_mat, n, ...
                                                r, Re, GM, normalized);
        [~, dUi, ~] = potentialGradient_nm(C_mat, S_mat, n - 1, ...
                                                r, Re, GM, normalized);

        Un(n, k) = vecnorm(dUf - dUi);
    end
end

% normalice values in a percentage
% % Un(1, :) = Un(1, :)./Un(1, 1).*100;
% % Un(2, :) = Un(2, :)./Un(2, 1).*100;
% % Un(3, :) = Un(3, :)./Un(3, 1).*100;
% % Un(4, :) = Un(4, :)./Un(4, 1).*100;
% % Un(5, :) = Un(5, :)./Un(5, 1).*100;
% % Un(6, :) = Un(6, :)./Un(6, 1).*100;

H = H./1000;


% plot 
pt = 2.5;
figure();
semilogy(H, Un(1, :), 'LineWidth', pt);
hold all;
semilogy(H, Un(2, :), 'LineWidth', pt);
semilogy(H, Un(3, :), 'LineWidth', pt);
semilogy(H, Un(4, :), 'LineWidth', pt);
semilogy(H, Un(5, :), 'LineWidth', pt);
semilogy(H, Un(6, :), 'LineWidth', pt);
xlabel('Orbit radius [km]');
ylabel('[m^2/s^2]');
grid on;
title('Gravity potential for different orbit radius')
legend('GM','n=2', 'n=3','n=4', 'n=5', 'n=6');
