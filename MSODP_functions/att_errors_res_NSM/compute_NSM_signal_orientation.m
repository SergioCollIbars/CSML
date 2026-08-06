clear;
clc;
close all;

addpath('functions/');
set(0,'defaultAxesFontSize',16);

%% USER SETTINGS

% Number of complete revolutions to plot
nRevolutions = 2;

% Mean Earth radius [km]
R_E = 6371.0088;

% Gravity-gradient data folder
folderPath = ...
    "/Users/sergiocollibars/Documents/GG_observations/120by120/500km/";

% Orbit file
orbitFile = "postion_velocity_A.mat";

%% MEASUREMENT MASK
%
%       xx xy xz yx yy yz zx zy zz

mask = logical([1, 0, 1, 0, 1, 0, 0, 0, 1]');

nSelected = sum(mask);

allComponentNames = [ ...
    "AA","AC","AR", ...
    "CA","CC","CR", ...
    "RA","RC","RR"];

selectedNames = allComponentNames(mask);

%% READ GRAVITY-GRADIENT OBSERVATIONS

GG_obs = parser_GG_obs_MSODP(folderPath);

referenceEpoch = datetime( ...
    2000,1,1,12,0,0, ...
    'TimeZone','UTC');

t_dateTime = referenceEpoch + seconds(GG_obs(:,1));

Nt = height(GG_obs);

%% ROTATE GRAVITY-GRADIENT OBSERVATIONS

[~,GG_nom] = rotate_GG_obs( ...
    GG_obs, ...
    zeros(3,Nt));

%% READ ORBIT DATA

orbitData = load(orbitFile).data;

orbitTime = datetime( ...
    orbitData(:,1), ...
    'ConvertFrom','juliandate', ...
    'TimeZone','UTC');

rOrbit = orbitData(:,2:4);
vOrbit = orbitData(:,5:7);

%% INTERPOLATE ORBIT AT THE GG EPOCHS

orbitJD = juliandate(orbitTime);
obsJD   = juliandate(t_dateTime);

% Sort orbit data
[orbitJD,sortIndex] = sort(orbitJD);

rOrbit = rOrbit(sortIndex,:);
vOrbit = vOrbit(sortIndex,:);

% Remove duplicate orbit epochs
[orbitJD,uniqueIndex] = unique(orbitJD,'stable');

rOrbit = rOrbit(uniqueIndex,:);
vOrbit = vOrbit(uniqueIndex,:);

% Find observation epochs covered by the orbit data
validTime = ...
    obsJD >= orbitJD(1) & ...
    obsJD <= orbitJD(end);

if ~all(validTime)

    warning( ...
        '%d observation epochs are outside the orbit-data interval.', ...
        sum(~validTime));

end

rECI = nan(Nt,3);
vECI = nan(Nt,3);

rECI(validTime,:) = interp1( ...
    orbitJD, ...
    rOrbit, ...
    obsJD(validTime), ...
    'linear');

vECI(validTime,:) = interp1( ...
    orbitJD, ...
    vOrbit, ...
    obsJD(validTime), ...
    'linear');

%% COMPUTE ORBIT LATITUDE

orbitLatitude = computeOrbitLatitudeECI(rECI,vECI);

% Convert to continuous angle
validOrbit = isfinite(orbitLatitude);

orbitLatitudeContinuousDeg = nan(Nt,1);

orbitLatitudeContinuousDeg(validOrbit) = rad2deg( ...
    unwrap(orbitLatitude(validOrbit)));

%% SELECT REQUESTED NUMBER OF REVOLUTIONS

firstValid = find(validOrbit,1,'first');

phaseStart = orbitLatitudeContinuousDeg(firstValid);
phaseEnd   = phaseStart + 360*nRevolutions;

plotMask = ...
    validOrbit & ...
    orbitLatitudeContinuousDeg >= phaseStart & ...
    orbitLatitudeContinuousDeg <= phaseEnd;

if ~any(plotMask)
    error('No orbit data are available in the requested interval.');
end

availablePhaseEnd = max(orbitLatitudeContinuousDeg(validOrbit));

if phaseEnd > availablePhaseEnd

    warning( ...
        ['Only %.2f revolutions are available. Plotting all ' ...
         'available data.'], ...
        (availablePhaseEnd-phaseStart)/360);

    phaseEnd = availablePhaseEnd;

    plotMask = ...
        validOrbit & ...
        orbitLatitudeContinuousDeg >= phaseStart & ...
        orbitLatitudeContinuousDeg <= phaseEnd;

end

fprintf( ...
    'Plotting %.2f orbital revolutions.\n', ...
    (phaseEnd-phaseStart)/360);

%% TIME FROM BEGINNING OF SELECTED INTERVAL

tHours = hours( ...
    t_dateTime - t_dateTime(firstValid));

%% COMPUTE EARTH-FIXED GROUND TRACK

rECEF = nan(Nt,3);

for j = find(validOrbit)'

    rECEF(j,:) = eciToEcefGMST( ...
        rECI(j,:), ...
        obsJD(j));

end

longitudeDeg = rad2deg(atan2( ...
    rECEF(:,2), ...
    rECEF(:,1)));

latitudeDeg = rad2deg(atan2( ...
    rECEF(:,3), ...
    hypot(rECEF(:,1),rECEF(:,2))));

longitudeDeg = mod(longitudeDeg + 180,360) - 180;

%% COMPUTE NSM WEIGHTS

disp('Computing NSM weights...');

weightContribution = nan( ...
    nSelected, ...
    Nt, ...
    nSelected);

for j = 1:Nt

    Hrot = compute_rotPartials_analy( ...
        GG_nom(j,2:end)', ...
        eye(3));

    % Apply measurement mask
    HrotSelected = Hrot(mask,:);

    % NSM basis vectors
    [~,~,V] = svd(HrotSelected');

    % Squared normalized components of each basis vector
    for k = 1:nSelected

        w = V(:,k);

        weightContribution(:,j,k) = ...
            abs(w).^2/sum(abs(w).^2);

    end

end

%% PLOT 1: EARTH-FIXED GROUND TRACK

figure;

plot( ...
    longitudeDeg(plotMask), ...
    latitudeDeg(plotMask), ...
    '.', ...
    'MarkerSize',8);

xlabel('Earth-fixed longitude [deg]');
ylabel('Geocentric latitude [deg]');

title(sprintf( ...
    'Ground track: %.2f revolutions', ...
    (phaseEnd-phaseStart)/360));

xlim([-180 180]);
ylim([-90 90]);

xticks(-180:60:180);
yticks(-90:30:90);

grid on;
box on;

%% PLOT 2: WEIGHT COMPONENTS VERSUS TIME

for k = 1:nSelected

    figure;
    hold on;

    for i = 1:nSelected

        weightCurrent = squeeze( ...
            weightContribution(i,:,k));

        plot( ...
            tHours(plotMask), ...
            weightCurrent(plotMask), ...
            'LineWidth',1.5, ...
            'DisplayName',selectedNames(i));

    end

    xlabel('Time from start [h]');
    ylabel('Weight contribution');

    title(sprintf( ...
        'NSM observable %d: weights versus time', ...
        k));

    ylim([0 1]);

    grid on;
    box on;

    legend( ...
        'Location','best', ...
        'NumColumns',2);

end

%% PLOT 3: WEIGHT COMPONENTS VERSUS ORBIT LATITUDE

% Orbit latitude wrapped to [0,360) degrees
orbitLatitudeDeg = mod( ...
    rad2deg(orbitLatitude), ...
    360);

for k = 1:nSelected

    figure;
    hold on;

    for i = 1:nSelected

        weightCurrent = squeeze( ...
            weightContribution(i,:,k));

        scatter( ...
            orbitLatitudeDeg(plotMask), ...
            weightCurrent(plotMask), ...
            18, ...
            'filled', ...
            'DisplayName',selectedNames(i));

    end

    xlabel('Orbit latitude [deg]');
    ylabel('Weight contribution');

    title(sprintf( ...
        'NSM observable %d: weights versus orbit latitude', ...
        k));

    xlim([0 360]);
    ylim([0 1]);

    xticks(0:45:360);

    grid on;
    box on;

    legend( ...
        'Location','best', ...
        'NumColumns',2);

end

%% PLOT 4: COMPONENT WEIGHTS OVER THE GROUND TRACK

for k = 1:nSelected

    figure;

    tiledlayout( ...
        2,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    for i = 1:nSelected

        nexttile;

        weightCurrent = squeeze( ...
            weightContribution(i,:,k));

        scatter( ...
            longitudeDeg(plotMask), ...
            latitudeDeg(plotMask), ...
            5, ...
            weightCurrent(plotMask), ...
            'filled');

        xlabel('Earth-fixed longitude [deg]');
        ylabel('Geocentric latitude [deg]');

        title(selectedNames(i));

        xlim([-180 180]);
        ylim([-90 90]);

        xticks(-180:60:180);
        yticks(-90:30:90);

        clim([0 1]);

        grid on;
        box on;

        cb = colorbar;
        cb.Label.String = 'Weight contribution';

    end

    sgtitle(sprintf( ...
        'NSM observable %d: component weights over ground track', ...
        k));

end
%% AUXILIARY FUNCTIONS

function orbitLatitude = computeOrbitLatitudeECI(rECI,vECI)
% Compute the argument of latitude in the instantaneous orbital plane.
%
% Output:
%   orbitLatitude: angle in the interval [0,2*pi)

    N = size(rECI,1);

    orbitLatitude = nan(N,1);

    kHat = [0;0;1];

    for j = 1:N

        r = rECI(j,:)';
        v = vECI(j,:)';

        if any(~isfinite(r)) || any(~isfinite(v))
            continue
        end

        h = cross(r,v);

        if norm(h) <= eps
            continue
        end

        hHat = h/norm(h);

        % Ascending-node direction
        nodeVector = cross(kHat,h);

        if norm(nodeVector) > 1e-12

            nodeHat = nodeVector/norm(nodeVector);

            % Second direction in the orbital plane
            qHat = cross(hHat,nodeHat);

            orbitLatitude(j) = atan2( ...
                dot(r,qHat), ...
                dot(r,nodeHat));

        else

            % Equatorial-orbit fallback
            orbitLatitude(j) = atan2(r(2),r(1));

        end

    end

    orbitLatitude = mod(orbitLatitude,2*pi);

end

function rECEF = eciToEcefGMST(rECI,julianDate)
% Approximate ECI-to-ECEF transformation using GMST.

    T = (julianDate - 2451545.0)/36525;

    gmstDeg = ...
        280.46061837 + ...
        360.98564736629*(julianDate - 2451545.0) + ...
        0.000387933*T.^2 - ...
        T.^3/38710000;

    theta = deg2rad(mod(gmstDeg,360));

    C_ECI2ECEF = [ ...
         cos(theta),  sin(theta), 0;
        -sin(theta),  cos(theta), 0;
         0,           0,          1];

    rECEF = (C_ECI2ECEF*rECI(:))';

end