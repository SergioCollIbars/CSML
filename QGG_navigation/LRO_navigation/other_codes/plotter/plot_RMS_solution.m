clear;
clc;
close all;
set(0, 'defaultAxesFontSize', 16);

addpath(['/Users/sergiocollibars/Desktop/CSML/QGG_navigation/' ...
         'LRO_navigation/other_codes/error_budget/']);

cspice_furnsh( ...
    '/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');

% Colors
errorColor = [0.5 0.5 0.5]; 
boundColor = 'b'; 

%% USER INPUTS

% Orbital period [s]
T_orbit = 9000;

% Monte Carlo realization to analyze
MC = 1;

%% INPUT FILES

folder = ...
    "/Users/sergiocollibars/Desktop/Lunar_orbit_simulator/data/";

files = [ ...
    "orbit_LRO_QGG_TS_10s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_20s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_30s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_40s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_50s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_60s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_120s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_180s_dataOut.mat", ...
    "orbit_LRO_QGG_TS_300s_dataOut.mat"];

Nfiles = numel(files);

%% PREALLOCATE RMS RESULTS

% Rows: files
% Columns: spacecraft-frame components [x, y, z]
rmsPosition = nan(Nfiles, 3);
rmsVelocity = nan(Nfiles, 3);

% Combined three-dimensional RMS
rmsPosition3D = nan(Nfiles, 1);
rmsVelocity3D = nan(Nfiles, 1);

% Optional RMS values of the 1-sigma covariance bounds
rmsSigmaPosition = nan(Nfiles, 3);
rmsSigmaVelocity = nan(Nfiles, 3);

%% PROCESS FILES

disp('Computing RMS values over the last orbital period ...');

for fl = 1:Nfiles

    filePath = fullfile(folder, files(fl));

    if ~isfile(filePath)
        warning('File not found: %s', filePath);
        continue;
    end

    loadedData = load(filePath);

    % Extract variables explicitly from the loaded file
    mC_struct  = loadedData.mC_struct;
    time       = loadedData.time;
    NB_MOON_mat = loadedData.NB_MOON_mat;

    % Ensure time is a row vector
    time = time(:)';

    fields = fieldnames(mC_struct);

    if MC > numel(fields)
        warning(['Requested MC realization %d does not exist in %s. ' ...
                 'Only %d realizations were found.'], ...
                 MC, files(fl), numel(fields));
        continue;
    end

    data = mC_struct.(fields{MC});

    Nt = numel(time);

    % Check consistency between time and state histories
    NtData = size(data.Xf, 2);

    if NtData ~= Nt
        warning(['Time and state lengths differ in %s. ' ...
                 'Using the minimum available length.'], files(fl));

        Nt = min(Nt, NtData);
        time = time(1:Nt);

        data.Xf         = data.Xf(:, 1:Nt);
        data.state_true = data.state_true(:, 1:Nt);
        data.Pf         = data.Pf(1:Nt, :);

        NB_MOON_mat = NB_MOON_mat(1:(3 * Nt), :);
    end

    % Determine samples in the last orbital period
    timeEnd   = time(end);
    timeStart = timeEnd - T_orbit;

    lastOrbitIndex = time >= timeStart & time <= timeEnd;

    if nnz(lastOrbitIndex) < 2
        warning(['Fewer than two samples were found in the final orbital ' ...
                 'period for %s.'], files(fl));
        continue;
    end

    actualWindow = time(end) - time(find(lastOrbitIndex, 1, 'first'));

    % Compute spacecraft orientation
    BN_mat = compute_orientation_SC( ...
        time, ...
        data.Xf(1:6, :)', ...
        data.I_ALIG, ...
        NB_MOON_mat);

    % Compute state errors and covariance
    [errorP, errorV, sigmaP, sigmaV] = ...
        get_pos_vel_cov(data, Nt, BN_mat);

    % RMS over only the final orbital period
    errorPLast = errorP(:, lastOrbitIndex);
    errorVLast = errorV(:, lastOrbitIndex);

    sigmaPLast = sigmaP(:, lastOrbitIndex);
    sigmaVLast = sigmaV(:, lastOrbitIndex);

    % Component-wise RMS:
    % rmsPosition(fl,:) = [RMS_x, RMS_y, RMS_z]
    rmsPosition(fl, :) = ...
        sqrt(mean(errorPLast.^2, 2, 'omitnan'))';

    rmsVelocity(fl, :) = ...
        sqrt(mean(errorVLast.^2, 2, 'omitnan'))';

    % Combined 3D RMS:
    %
    % sqrt(mean(ex^2 + ey^2 + ez^2))
    %
    % This is also equal to:
    % sqrt(RMS_x^2 + RMS_y^2 + RMS_z^2)
    rmsPosition3D(fl) = ...
        sqrt(mean(sum(errorPLast.^2, 1), 'omitnan'));

    rmsVelocity3D(fl) = ...
        sqrt(mean(sum(errorVLast.^2, 1), 'omitnan'));

    % Optional: RMS of the component-wise 3-sigma values
    rmsSigmaPosition(fl, :) = ...
        3.*sqrt(mean(sigmaPLast.^2, 2, 'omitnan'))';

    rmsSigmaVelocity(fl, :) = ...
        3.*sqrt(mean(sigmaVLast.^2, 2, 'omitnan'))';

    %% Display results

    fprintf('\nFile %d: %s\n', fl, files(fl));
    fprintf('Final-orbit window: %.2f s\n', actualWindow);
    fprintf('Number of samples: %d\n', nnz(lastOrbitIndex));

    fprintf('Position RMS components:\n');
    fprintf('  X: %.6e\n', rmsPosition(fl, 1));
    fprintf('  Y: %.6e\n', rmsPosition(fl, 2));
    fprintf('  Z: %.6e\n', rmsPosition(fl, 3));
    fprintf('  Combined 3D RMS: %.6e\n', rmsPosition3D(fl));

    fprintf('Velocity RMS components:\n');
    fprintf('  X: %.6e\n', rmsVelocity(fl, 1));
    fprintf('  Y: %.6e\n', rmsVelocity(fl, 2));
    fprintf('  Z: %.6e\n', rmsVelocity(fl, 3));
    fprintf('  Combined 3D RMS: %.6e\n', rmsVelocity3D(fl));
end

%% STORE RESULTS IN A TABLE

resultsTable = table( ...
    files(:), ...
    rmsPosition(:, 1), ...
    rmsPosition(:, 2), ...
    rmsPosition(:, 3), ...
    rmsPosition3D, ...
    rmsVelocity(:, 1), ...
    rmsVelocity(:, 2), ...
    rmsVelocity(:, 3), ...
    rmsVelocity3D, ...
    'VariableNames', { ...
        'File', ...
        'PositionRMS_X', ...
        'PositionRMS_Y', ...
        'PositionRMS_Z', ...
        'PositionRMS_3D', ...
        'VelocityRMS_X', ...
        'VelocityRMS_Y', ...
        'VelocityRMS_Z', ...
        'VelocityRMS_3D'});

disp(' ');
disp('RMS results for the last orbital period:');
disp(resultsTable);

%% PLOT RESULTS
figure();
fs = 1/60; R = 1/6;
% % xAxis = [5E-4 /fs / R , 2.7E-3 /fs / R, 8E-3 / fs / R];
xAxis = [10, 20, 30, 40, 50, 60, 120, 180, 300];
axis  = ["Radial", "Tangential", "Normal"];
for k = 1:3
    subplot(3, 1, k)
    plot(xAxis, rmsSigmaPosition(:, k)', ...
        'Marker','square', 'MarkerFaceColor', 'b');
    title(axis(k))
     grid on; xlabel('sampling time [sec]'); ylabel('[m]');
     % % ylim([0, 20]);
end


%% CLEAR SPICE KERNELS

cspice_kclear;


%% AUXILIARY FUNCTIONS

function [errorP, errorV, sigmaP, sigmaV] = ...
    get_pos_vel_cov(data, Nt, BN_mat)

    % State error in the inertial frame
    error = data.state_true - data.Xf;

    sigmaP = nan(3, Nt);
    sigmaV = nan(3, Nt);
    errorP = nan(3, Nt);
    errorV = nan(3, Nt);

    for k = 1:Nt

        % Reconstruct covariance matrix
        p = reshape(data.Pf(k, :), [18, 18]);

        maxInd = 3 * k;
        minInd = maxInd - 2;

        % Inertial-to-spacecraft-frame rotation
        BN = BN_mat(minInd:maxInd, :);

        % Transform position and velocity covariance blocks
        A = blkdiag(BN, BN, eye(6), eye(6));

        P_B = A * p * A';

        % Ensure numerical roundoff does not produce a small negative
        % value inside the square root
        variance = max(diag(P_B), 0);
        stateStd = sqrt(variance);

        sigmaP(:, k) = stateStd(1:3);
        sigmaV(:, k) = stateStd(4:6);

        % Transform position and velocity errors
        errorP(:, k) = BN * error(1:3, k);
        errorV(:, k) = BN * error(4:6, k);
    end
end