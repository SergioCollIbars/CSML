clear;
clc;
close all;
delete(findall(groot,'Type','figure'))

addpath("data/")
addpath(genpath("functions/"))

set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm');
%%                      LOW LUNAR NAVIGATION CODE
%
% Date: 15/12/2025
% Author: Sergio Coll-Ibars
% Description: Using simulated low lunar trajectories and gradiometer
% measurements, estimate the position and velocity of the S/C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load meta data
metaData_path = "Metadata.txt";
mtd = readParams("data/"+metaData_path);

% start Simulation
[planetParams, Cnm_list, Snm_list, ...
    state0, t_range]                    = loadUniverse(metaData_path);
[instrumentParams, instrument_alig]     = loadInstrument(metaData_path);

% integrate Trajectory
disp('Simulating trajectory ...')
tmin = t_range(1);
tmax = t_range(2);
t = linspace(tmin, tmax, (tmax-tmin) * instrumentParams(1, 5));

options = odeset('RelTol',1e-13,'AbsTol',1e-13); Nx = 6;
PHI0    = reshape(eye(Nx,Nx), [Nx*Nx, 1]);
X0      = state0;
[time, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
    Cnm_list, Snm_list), t, [X0;PHI0], ...
    options);

disp('  DONE ...')

% compute S/C orientation
disp('Computing planet orientation ...')
[NB_EARTH_mat, NB_MOON_mat] = compute_orientation_planets(time);
disp('  DONE ...')

% generate Measurements (Inertial frame)
disp('Simulating measurements ...')
[Y, bias, signal_err] = compute_measurements(instrumentParams, planetParams, ...
    time, state, Cnm_list, Snm_list, NB_EARTH_mat, NB_MOON_mat);

disp('  DONE ...')

% Filtering process
disp('Running Filter ...')
state_true = [state(:, 1:6)';bias];
[Xf, Pf]   = filter_measurements(metaData_path,time,state_true,...
    instrument_alig, NB_EARTH_mat, NB_MOON_mat, Cnm_list, Snm_list, Y, ...
    signal_err, planetParams, instrumentParams);
disp('  DONE ...')

% Plot results
if(mtd.plot)
    disp('Plotting Results ...')
    plot_trajectory(time, state);
    
    plot_measurements(time, Y, bias, signal_err, instrumentParams,...
        Xf, instrument_alig); 
    
    mask = instrumentParams(:, 1);
    plot_results(time, state, bias, Xf, Pf, mask, instrument_alig);
    disp('  DONE ...')
end

% Save results
if(mtd.save)
    disp('Saving Results ...')
    save_results(Xf, Pf, mtd.folder);
    disp('  DONE ...')
end

% clear kernels
cspice_kclear
