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

for fld = 1:length(mtd.folder)
    folder_name = string(mtd.folder{fld});
    disp('Running ' + folder_name + '(' +...
        string(fld) + '/' + string(length(mtd.folder)) + ')');

    % start Simulation
    [planetParams, Cnm_list, Snm_list, ...
        state0, t_range]              = loadUniverse(mtd.folder{fld});
    [instrumentParams_GG]             = loadInstrument_GG(mtd.folder{fld});
    [instrumentParams_ST]             = loadInstrument_ST(mtd.folder{fld});
    
    % integrate Trajectory
    disp('Simulating trajectory ...')
    tmin = t_range(1);
    tmax = t_range(2);
    t = linspace(tmin, tmax, (tmax-tmin) * instrumentParams_GG(1, 5));
    
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
    [Y, bias, signal_err, SF] = compute_measurements(instrumentParams_GG, ...
        planetParams, time, state, Cnm_list, Snm_list,...
        NB_EARTH_mat, NB_MOON_mat);
    
    disp('  DONE ...')
    
    % Filtering process
    disp('Running Filter ...')
    state_true = [state(:, 1:6)';bias; SF];
    [Xf, Pf, posfit, I_ALIG]   = filter_measurements(mtd.folder{fld},time,state_true,...
        NB_EARTH_mat, NB_MOON_mat, Cnm_list, Snm_list, Y, ...
        signal_err, planetParams, instrumentParams_GG, instrumentParams_ST);
    disp('  DONE ...')
    
    % Plot results
    if(mtd.plot)
        disp('Plotting Results ...')
        plot_trajectory(time, state, folder_name);
        
        plot_measurements(time, Y, bias, signal_err, instrumentParams_GG,...
            Xf, I_ALIG, folder_name, NB_MOON_mat, SF); 
        
        mask = instrumentParams_GG(:, 1);
        plot_results(time, state, bias, SF, Xf, Pf, posfit, mask, ...
            I_ALIG, folder_name, NB_MOON_mat);
        disp('  DONE ...')
    end
    
    % Save results
    if(mtd.save)
        disp('Saving Results ...')
        save_results(Xf, Pf, mtd.folder{fld});
        disp('  DONE ...')
    end
end

% clear kernels
cspice_kclear
