clear;
clc;
close all;
format long g;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
%%                    QGG NAVIGATION OBSERVABILITY
% Description: compute the information matrix and upper bpund at different
% locations in space.
% Author: Sergio Coll Ibars
% Date: 11/12/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%                COMPUTE INFO AT DIFFERENT LONG AND LAT

% define system
tmin = 0;                                          % Phi = 0  [-]
tmax = 3.4968 + tmin;                              % [-]
frec = 1/60;

system = "CR3BP"; % options: 2BP, CR3BP, F2BP
[planetParams, poleParams, C_mat, S_mat, ~, ~] = ...
    load_universe(system, [tmin, tmax], frec);

% Domain definiton. Non-dimensional units 
Re = 6371E3;                    % Earth radii[m]
Re_ND = Re / planetParams(2);   % [-]
Rm = 1740E3;                    % Moon radii [m]
Rm_ND = Rm / planetParams(2);   % [-]
MoonD  = 1 - planetParams(1);   % [-]
EarthD = - planetParams(1);     % [-]

% spherical coordinates
N = 100;
r = 2*Re_ND;                            % [-]
PHI = linspace(-pi/2, pi/2, N);         % [rad]
LAMBDA = linspace(0, pi, N);            % [rad]

BD_QGG = ones(N, N);
for j = 1:N     % phi values
    for k = 1:N % lambda values
        phi = PHI(j);
        lambda = LAMBDA(k);
        BD_QGG(k, j) = compute_BD(r, phi, lambda, planetParams, poleParams, C_mat, S_mat, system);
    end
end

figure();
[X, Y] = meshgrid(rad2deg(PHI), rad2deg(LAMBDA));
surf(X, Y, real(BD_QGG), 'EdgeColor','none')
xlabel('Latitude \phi [deg]')
ylabel('Longitude \lambda [deg]')
h = colorbar;
h.Ticks = linspace(min(caxis), max(caxis), 5); % Set tick positions
h.TickLabels = {'Low', 'Medium-Low', 'Medium', 'Medium-High', 'High'}; % Custom labels
xlim([-90, 90])


%%                  COMPUTE INFO. AT DIFFERENT POINTS


% meshgrid parameters
Ref= Re_ND;
N = 150;
% % valsX = [-1.5*Ref + EarthD, 1.5*Ref + EarthD, N];
valsX = [-150*Ref, 150*Ref, N];
valsY = [-150*Ref, 150*Ref, N];
valsZ = [-150*Ref, 150*Ref, N];

[IF_QGG, BD_QGG, X, Y, Z] = compute_boundingContour(valsX, valsY, valsZ,...
    planetParams, poleParams, C_mat, S_mat, system, Ref);

% projection in the Y = 0 plane
[~, index] = min(abs(Y(:, 1, 1) - 0));
X_slice = squeeze(X(index, :, :));
Z_slice = squeeze(Z(index, :, :));
IF_QGG_slice = squeeze(IF_QGG(index, :, :));
BD_QGG_slice = squeeze(BD_QGG(index, :, :));

figure;
contourf(X_slice./Ref, Z_slice./Ref, log10(BD_QGG_slice), 20, 'LineWidth', 1.2);
xlabel('X');
ylabel('Z');
colormap('turbo')
title('Contour Projection on Y = 0 Plane');
axis equal
colorbar

% projection in the Z = 0 plane
[~, index] = min(abs(Z(1, 1, :) - 0));
Y_slice = squeeze(Y(:, :, index));
X_slice = squeeze(X(:, :, index));
IF_QGG_slice = squeeze(IF_QGG(:, :, index));
BD_QGG_slice = squeeze(BD_QGG(:, :, index));

figure;
contourf(X_slice./Ref, Y_slice./Ref, log10(BD_QGG_slice), 20, 'LineWidth', 1.2);
xlabel('X');
ylabel('Y');
colormap('turbo')
title('Contour Projection on Z = 0 Plane');
axis equal
colorbar

scale  = planetParams(2)/1E3;   % [km]
lg = ["", "",  "Earth", "Moon", "L_1", "L_2", "L_3", "L_4", "L_5"];
ttitle = "Information gain value (\delta |\Delta_x|) in log10 scale";
plot_sensitivityContour(X, Y, log10(IF_QGG), ttitle, ...
    planetParams(1), scale, lg);

ttitle = "PCRLB upper bound around Moon. Gradiometer measurement";
plot_sensitivityContour(X, Y, log10(BD_QGG), ttitle, ...
    planetParams(1), scale, lg);


% clear kernels
cspice_kclear

%% FUNCTIONS
function [] = plot_sensitivityContour(X, Y, Z, ttitle, mu, scale, lg)
    % Lagrange points
    L_1 = .836915 * scale;
    L_2 = 1.15568 * scale;
    L_3 = -1.00506 * scale;
    
    % trajectory
    figure()
    contourf(X.*scale, Y.*scale, Z, 'EdgeColor', 'none');
    view(0, 90);
    colorbar()
    xlabel('X [Km]')
    ylabel('Y [Km]')
    hold all;
    axis equal;
    plot(-mu*scale,0, "o",'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
    plot((1-mu)*scale,0,"pentagram",'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
    plot(L_1,0,'rv','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot(L_2,0,'r^','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot(L_3,0,'rp','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot((0.5-mu)*scale,(sqrt(3)/2)*scale,'rX','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot((0.5-mu)*scale,(-sqrt(3)/2)*scale,'rs','MarkerFaceColor','k', 'MarkerEdgeColor','k')
    legend('Earth', 'Moon', 'L1', 'L2', 'L3', 'L4', 'L5')
    legend(lg)
    title(ttitle)
end

function [IF_QGG, BD_QGG, X, Y, Z] = compute_boundingContour(valsX, ...
    valsY, valsZ, planetParams, poleParams, C_mat, S_mat, system, Ref)
    
    % Define output matrices
    IF_QGG   = ones(valsX(3), valsY(3), valsZ(3)) * NaN;
    BD_QGG   = ones(valsX(3), valsY(3), valsZ(3)) * NaN;
    Nx = valsX(3); Ny = valsY(3); Nz = valsZ(3);
    
    % Define Contour domain
    [X, Y, Z] = meshgrid(linspace(valsX(1), valsX(2), valsX(3)), linspace(valsY(1), valsY(2), valsY(3)), ...
        linspace(valsZ(1), valsZ(2), valsZ(3)));
    
    % measurements matrices
    sigmaMeas = [1, 1/2] * 1E-12;                           % [1/s^2]
    sigmaMeas = sigmaMeas./(planetParams(3)^2);             % [-]
    R_QGG = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
        sigmaMeas(2), sigmaMeas(1)].^2);                    % [-]
    
    for j = 1:Nx    % x domain
        for i =1:Ny % y domain
            for k = 1:Nz
                time = 0;   % [-]
                state = [X(i, j, k), Y(i, j, k), Z(i, j, k)];
        
                % gradiometer sensitivity
                [posE, posM] = compute_posPrimaries(time, planetParams, system);
        
                [~, Hc, ~] = compute_measurements(time, state, planetParams, poleParams, ...
                    C_mat, S_mat, 0, 0, 0, [], 0, posE, posM, [], system);
        
                % visibility matrix
                h = Hc(:, 1:3);

                relPos_E = vecnorm(state' - posE);
                relPos_M = vecnorm(state' - posM);
        
                % rescale measurements
                scale  = planetParams(2);   % [m]
        
                f = h' * inv(R_QGG) * h;
                l = min(eig(f));
                BD = ((1/l))^(1/2);                              % bounding for largest uncertanty direction
                IF  = det(f);

                if(relPos_M < Ref || relPos_E < Ref)
                    BD_QGG(i, j, k)   = NaN;                  % [m]
                    IF_QGG(i, j, k)   = NaN;                  % [1/m^3]
                else
                    BD_QGG(i, j, k)   = BD * scale;           % [m]
                    IF_QGG(i, j, k)   = IF * 1/scale^3;       % [1/m^3]
                end
            end
        end
    end
end

function [BD_QGG] = compute_BD(r, phi, lambda, planetParams, poleParams, C_mat, S_mat, system)
    % measurements matrices
    sigmaMeas = [1, 1/2] * 1E-12;                           % [1/s^2]
    sigmaMeas = sigmaMeas./(planetParams(3)^2);             % [-]
    R_QGG = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
        sigmaMeas(2), sigmaMeas(1)].^2);                    % [-]

    time = 0;   % [-]
    x = r * cos(phi) * sin(lambda);
    y = r * cos(phi) * cos(lambda);
    z = r * sin(phi);
    pos = [x, y, z];

    % gradiometer sensitivity
    [posE, posM] = compute_posPrimaries(time, planetParams, system);

    state = pos + posE';

    [~, Hc, ~] = compute_measurements(time, state, planetParams, poleParams, ...
        C_mat, S_mat, 0, 0, 0, [], 0, posE, posM, [], system);

    % visibility matrix
    h = Hc(:, 1:3);

    % rescale measurements
    scale  = planetParams(2);   % [m]

    f = h' * inv(R_QGG) * h;
    l = min(eig(f));
    BD = ((1/l))^(1/2);                              % bounding for largest uncertanty direction

    BD_QGG   = BD * scale;
end