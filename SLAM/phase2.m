clear;
clc;
close all;

addpath("functions/"); addpath("data/");

set(0,'defaultAxesFontSize',16);

%%                       PHASE 2 ANALYSIS
% Description: Simulate SLAM process during the Phase 0 of the mission.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load meta data
metaData_path = "Metadata.txt";
mtd = readParams("data/"+metaData_path);

% start Simulation
[planetParams, poleParams, Cnm, Snm, ...
    initCond, t_range]                 = loadUniverse(metaData_path);
[instrumentParams]                     = loadInstrument(metaData_path);

% integrate Trajectory
disp('Simulating trajectory ...')
tmin = t_range(1);
tmax = t_range(2);
t = linspace(tmin, tmax, tmax * instrumentParams(1, 5));

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, ~, Ncs] = count_num_coeff(planetParams(3));
Nx      = 12 + Ncs; 
PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
[time, state] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
    Cnm, Snm), t, [initCond(1:3);initCond(4:6);...
    instrumentParams(:, 3);PHI0], options);

plot_trajectory(time, state)
disp('  DONE ...')

% compute S/C orientation
disp('Computing S/C orientation ...')
[BN_mat] = compute_orientation(time, state, instrumentParams);
disp('  DONE ...')

% generate Measurements
disp('Simulating measurements ...')
[Y, state] = compute_measurements(instrumentParams, planetParams, ...
    poleParams, time, state, Cnm, Snm, BN_mat);

plot_measurements(time, Y)
disp('  DONE ...')

% SLAM process
disp('Running SLAM ...')
[Xs, Ps, Xg, Pg, h, ax] = do_SLAM_v3(metaData_path, time, state, ...
    planetParams, poleParams, instrumentParams, BN_mat, Cnm, Snm, Y);
disp('  DONE ...')

% % % Smoothing
% % [Xg, Pg] = do_smoothing(time, state, planetParams, poleParams,...
% %     instrumentParams, BN_mat, Cnm, Snm, Y, Xs, Ps, Xg, Pg, h, ax);
% % disp('  DONE ...')

% save gravity field data
if(mtd.save)
    disp('Save Gravity data ...');
    save("data/"+mtd.folder+"/gravField_cov.mat", 'Pg');
    save("data/"+mtd.folder+"/gravField_coeff.mat", 'Xg');
end