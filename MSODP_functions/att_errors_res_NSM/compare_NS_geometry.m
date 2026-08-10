clear;
clc;
close all;

%% COMPARE NULL SPACE GEOMETRY
% Description: Using the GG observations generated in MSODP, compute the
% null space geometry in measurement domain.
%
% Author: Sergio Coll-Ibars
% Date: 08/10/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Measurement Mask
%               xx xy xz yx yy yz zx zy zz
mask = logical([1, 0, 1, 0, 1, 0, 0, 0, 1]');

allComponentNames = [ ...
    "AA","AC","AR", ...
    "CA","CC","CR", ...
    "RA","RC","RR"];

componentNames = allComponentNames(mask);


%% Extract GG observations

% 120x120 gravity field
folderPath = ...
    "/Users/sergiocollibars/Documents/GG_observations/120by120/500km";

GG = parser_GG_obs_MSODP(folderPath);

[~, GG_nom_120x120] = rotate_GG_obs( ...
    GG, zeros(3, length(GG(:,1))));


% 5x5 gravity field
folderPath = ...
    "/Users/sergiocollibars/Documents/GG_observations/5by5/";

GG = parser_GG_obs_MSODP(folderPath);

[~, GG_nom_5x5] = rotate_GG_obs( ...
    GG, zeros(3, length(GG(:,1))));


clear GG folderPath;

time = GG_nom_5x5(:,1) - GG_nom_5x5(1,1);


%% Compute NS geometry

disp('Compute Null Space geometry ...')

nMeas = sum(mask);
Nt    = length(time);

w_0x0     = nan(nMeas, nMeas, Nt);
w_5x5     = nan(nMeas, nMeas, Nt);
w_120x120 = nan(nMeas, nMeas, Nt);

% Attitude-error partial matrix for a Point Mass gravity field
% Rows:    [AA AC AR CA CC CR RA RC RR]
% Columns: [dtheta_A dtheta_C dtheta_R]
Hrot_0x0 = 3* [ ...
    0,  0,  0;   % AA
    0,  0,  0;   % AC
    0,  1,  0;   % AR
    0,  0,  0;   % CA
    0,  0,  0;   % CC
    -1,  0,  0;  % CR
    0,  1,  0;   % RA
    -1,  0,  0;   % RC
    0,  0,  0];  % RR

for j = 1:Nt

    % GG attitude error partials
    Hrot_5x5 = compute_rotPartials_analy( ...
        GG_nom_5x5(j,2:end)', ...
        eye(3));

    Hrot_120x120 = compute_rotPartials_analy( ...
        GG_nom_120x120(j,2:end)', ...
        eye(3));

    % Apply measurement mask
    HrotSelected_0x0 = Hrot_0x0(mask,:);

    HrotSelected_5x5 = Hrot_5x5(mask,:);

    HrotSelected_120x120 = Hrot_120x120(mask,:);

    % NSM basis vectors
    [~,~,V_0x0] = svd(HrotSelected_0x0');

    [~,~,V_5x5] = svd(HrotSelected_5x5');

    [~,~,V_120x120] = svd(HrotSelected_120x120');

    % Squared normalized components of each basis vector
    w_0x0(:,:,j)     = V_0x0.^2;
    w_5x5(:,:,j)     = V_5x5.^2;
    w_120x120(:,:,j) = V_120x120.^2;

end


%% Difference in NS geometry

w_diff_5_and_120 = w_5x5 - w_120x120;
w_diff_0_and_120 = w_0x0 - w_120x120;


%% Plot NS geometry: S2

plot_weights( ...
    time, ...
    w_0x0, ...
    2, ...
    componentNames, ...
    'S_2 for a 0 by 0 grav. field');

plot_weights( ...
    time, ...
    w_5x5, ...
    2, ...
    componentNames, ...
    'S_2 for a 5 by 5 grav. field');

plot_weights( ...
    time, ...
    w_120x120, ...
    2, ...
    componentNames, ...
    'S_2 for a 120 by 120 grav. field');

plot_weights( ...
    time, ...
    w_diff_5_and_120, ...
    2, ...
    componentNames, ...
    'S_2 weight difference: 5x5 - 120x120');

plot_weights( ...
    time, ...
    w_diff_0_and_120, ...
    2, ...
    componentNames, ...
    'S_2 weight difference: 0x0 - 120x120');


%% Plot NS geometry: S3

plot_weights( ...
    time, ...
    w_0x0, ...
    3, ...
    componentNames, ...
    'S_3 for a 0 by 0 grav. field');

plot_weights( ...
    time, ...
    w_5x5, ...
    3, ...
    componentNames, ...
    'S_3 for a 5 by 5 grav. field');

plot_weights( ...
    time, ...
    w_120x120, ...
    3, ...
    componentNames, ...
    'S_3 for a 120 by 120 grav. field');

plot_weights( ...
    time, ...
    w_diff_5_and_120, ...
    3, ...
    componentNames, ...
    'S_3 weight difference: 5x5 - 120x120');

plot_weights( ...
    time, ...
    w_diff_0_and_120, ...
    3, ...
    componentNames, ...
    'S_3 weight difference: 0x0 - 120x120');


%% AUXILIARY FUNCTIONS

function plot_weights(time, weights, basisIndex, componentNames, plotTitle)
% plot_weights
%
% Plot the contribution of each gravity-gradient measurement component
% to a selected NSM basis vector.
%
% Inputs:
%   time           - Time vector [s]
%   weights        - Weight matrix [nMeas x nBasis x Nt]
%   basisIndex     - Basis vector to plot (e.g., 2 for S2)
%   componentNames - Names of selected GG components
%   plotTitle      - Figure title

    data = squeeze(weights(:,basisIndex,:))';

    figure();
    hold on;

    for k = 1:size(data,2)

        plot( ...
            time./3600, ...
            data(:,k), ...
            'LineWidth',2, ...
            'DisplayName',componentNames(k));

    end

    grid on;
    xlim([0,3]);

    xlabel('Time [hours]');
    ylabel('Weight');

    title(plotTitle);

    legend( ...
        'show', ...
        'Location','best');

end