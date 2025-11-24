clear;
clc;
close all;

addpath("functions/"); addpath("data/");

set(0,'defaultAxesFontSize',16);

%%                       PHASE 0 ANALYSIS
% Description: Simulate SLAM process during the Phase 0 of the mission.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start Simulation
[planetParams, poleParams, Cnm, Snm, ...
    initCond]                          = loadUniverse();
[instrumentParams]                     = loadInstrument();

% integrate Trajectory
disp('Simulating trajectory ...')
tmin = 0;
tmax = 1 * 86400;
t = linspace(tmin, tmax, tmax * instrumentParams(1, 5));

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, ~, Ncs] = count_num_coeff(planetParams(3));
Nx      = 12 + Ncs; 
PHI0 = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
[time, state] = ode113(@(t, x) EoM(t, x, planetParams, poleParams, ...
    Cnm, Snm), t, [initCond(1:3);initCond(4:6);...
    instrumentParams(:, 3);PHI0], options);

plot_trajectory(time, state)

% compute S/C orientation
disp('Computing S/C orientation ...')
[BN_mat] = compute_orientation(time, state, instrumentParams);

% generate Measurements
disp('Simulating measurements ...')
[Y] = compute_measurements(instrumentParams, planetParams, ...
    poleParams, time, state, Cnm, Snm, BN_mat);

plot_measurements(time, Y)

% SLAM process
disp('Running SLAM ...')
[Xs, Ps, Xg, Pg] = do_SLAM(time, state, planetParams, ...
    poleParams, instrumentParams, BN_mat, Cnm, Snm, Y);

% Smoothing
do_smoothing(time, state, planetParams, poleParams, instrumentParams, ...
    BN_mat, Cnm, Snm, Y, Xs, Ps, Xg, Pg, h, ax);