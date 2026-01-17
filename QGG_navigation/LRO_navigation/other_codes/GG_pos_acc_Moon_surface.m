clear; clc; close all;
format long g;

addpath('../data/'); addpath(genpath("../functions/"));
set(0,'defaultAxesFontSize',16);
%%          COMPUTE POSITION UNCERTAINTY FROM A SINGLE MEASUREMENT
% Description: At a given altitude, compute the position uncertianty from a
% single measurement over different spatial locations
% Author: Sergio Coll-Ibars
% Date: 01/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gravity field
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";
input_coeffuncrt = "COEFSUNCRT_MOON_GRGM1200.txt";

file             = readmatrix(input_coeffuncrt);
R_M = file(2);   normalized = file(3); GM_moon = 4.9028001224453001e+03 * 1E9;
n_max            = 100; 
params           = [GM_moon, R_M, normalized, n_max];
[Nc, Ns, Ncs]    = count_num_coeff(n_max);

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);

file             = readmatrix(input_coeffuncrt);
SH_uncrt         = file(4:end);

[Cnm, Snm] = list2mat(n_max, Nc, Ns, SH_coeff);

%% Meshgrid

% Resolution
dlat = 0.3;      % deg
dlon = 0.3;      % deg
Ncountour = 50;

lat = (-90:dlat:90)';         % column [deg]
lon = (-180:dlon:180);        % row    [deg]

[Lon, Lat] = meshgrid(lon, lat);   % 2D grids (size: nLat x nLon)
[count, unct] = compute_Moon_contour(Lon,Lat, params, Cnm, Snm);

% plot
figure();
contourf(Lon, Lat, count, Ncountour, 'LineColor','none');
hold on
contour(Lon, Lat, count, Ncountour, 'LineWidth', 0.1, 'Color', 'k');

maxVal = max(max(max(log10(unct))));
minVal = min(min(min(log10(unct))));
tt     =["E", "N", "U"];
for k = 1:3
    figure()
    Z = log10(squeeze(unct(:, :, k)));
    
    contourf(Lon, Lat, Z, Ncountour, 'LineColor','none');
    colorbar; % Lock colormap limits
    colormap(turbo);
    clim([minVal maxVal]);
    
    hold on
    contour(Lon, Lat, count, Ncountour, 'LineWidth', 0.1, 'Color', 'k');
    xlabel('\lambda [deg]'); ylabel('\phi [deg]');
    title(tt(k));
end

figure()
Z = log10(squeeze(unct(:, :, 4)));

contourf(Lon, Lat, Z, Ncountour, 'LineColor','none');
colorbar; % Lock colormap limits
colormap(turbo);
clim([minVal maxVal]);

hold on
contour(Lon, Lat, count, Ncountour, 'LineWidth', 0.01, 'Color', 'k');
xlabel('\lambda [deg]'); ylabel('\phi [deg]');
title('Position uncertainty norm');

