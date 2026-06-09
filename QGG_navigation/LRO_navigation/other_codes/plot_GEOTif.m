clear; clc; close all;

%% Read GeoTIFF
filename = 'gggrx_1200a_anomerr_l660.tif';

color_limits = [0 5];   % mGal, as in the paper-style figure

%% ------------------------------------------------------------------------
% Read GeoTIFF
% -------------------------------------------------------------------------
[A, R] = readgeoraster(filename);
A = double(A);

% Remove common no-data values if present
A(A < -1e20) = NaN;
A(A >  1e20) = NaN;

% Optional: clip values for visualization
A(A < color_limits(1)) = color_limits(1);
A(A > color_limits(2)) = color_limits(2);

%% ------------------------------------------------------------------------
% Construct latitude/longitude using PDS storage convention
%
% Product convention:
%   Projection: simple cylindrical
%   Longitude: 0 to 360 deg East
%   Latitude:  -90 to 90 deg
%   Resolution: 16 pixels/degree
%
% For the 1200A map:
%   nlat = 2880
%   nlon = 5760
% -------------------------------------------------------------------------
nlat = size(A, 1);
nlon = size(A, 2);

dlat = 180 / nlat;
dlon = 360 / nlon;

% Pixel-center coordinates
lat = 90 - dlat/2 - (0:nlat-1) * dlat;    % north to south
lon = 270  + dlon/2 + (0:nlon-1) * dlon;    % 0 to 360 East

%% ------------------------------------------------------------------------
% Convert longitude from 0...360 to -180...180 and reorder raster columns
% -------------------------------------------------------------------------
lon_wrapped = lon;
lon_wrapped(lon_wrapped > 180) = lon_wrapped(lon_wrapped > 180) - 360;

[lon_plot, idx_lon] = sort(lon_wrapped);
A_plot = A(:, idx_lon);

lat_plot = lat;

[Lon, Lat] = meshgrid(lon_plot, lat_plot);

%% ------------------------------------------------------------------------
% Diagnostic plot in regular latitude/longitude axes
% -------------------------------------------------------------------------
figure('Color','w');

imagesc(lon_plot, lat_plot, A_plot);
set(gca, 'YDir', 'normal');

axis tight;
grid on;

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('Gravity Anomaly Error: Latitude/Longitude Check');

colormap(jet);
caxis(color_limits);

cb = colorbar;
ylabel(cb, 'mGal');

%% ------------------------------------------------------------------------
% Mollweide projection plot
% -------------------------------------------------------------------------
figure('Color','w');

axesm('mollweid', ...
    'MapLatLimit', [-90 90], ...
    'MapLonLimit', [-180 180], ...
    'Frame', 'on', ...
    'Grid', 'on', ...
    'MeridianLabel', 'off', ...
    'ParallelLabel', 'on', ...
    'PLabelLocation', [-60 -30 0 30 60], ...
    'MLineLocation', 30, ...
    'PLineLocation', 30, ...
    'FontSize', 14, ...
    'LabelFormat', 'none');

axis off;

surfm(Lat, Lon, A_plot);
shading interp;

colormap(jet);
caxis(color_limits);

tightmap;

% Grid and frame style
setm(gca, ...
    'GLineStyle', ':', ...
    'GLineWidth', 0.5, ...
    'GColor', [0.25 0.25 0.25], ...
    'FLineWidth', 1.2);

% Adjust map position
ax = gca;
ax.Position = [0.07 0.08 0.86 0.72];

% Colorbar like the paper
cb = colorbar('northoutside');
cb.Label.String = 'mGal';
cb.Label.FontSize = 14;
cb.FontSize = 12;
cb.Position = [0.14 0.83 0.72 0.035];

title('Gravity Anomaly Error', ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% Optional save
% -------------------------------------------------------------------------
% exportgraphics(gcf, 'GRAIL_gravity_anomaly_error_mollweide.png', ...
%     'Resolution', 300);