function out = run_one_MC_spice(metaKernelPath)
    cspice_furnsh(metaKernelPath);
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
    end

    out = Xf;
end
