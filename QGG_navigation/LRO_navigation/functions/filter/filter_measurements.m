function [X_EKF, P_EKF] = filter_measurements(metaData_file,time,state_true, ...
    instrument_alig, NB_EARTH_mat, NB_MOON_mat, Cnm_true, Snm_true, ...
    Y_true, signal_err, planetParams, instrumentParams)
    %% FILTER GRADIOMETER MEASUREMENTS
    
    % load filter parameters
    [R0, P0, Q0, Qb, delta_state0, Cnm_filt, Snm_filt] = ...
        loadFilterParams(metaData_file, planetParams, instrumentParams,...
        Cnm_true, Snm_true);

    % get measurement mask
    mask = instrumentParams(:, 1);
    fs   = instrumentParams(1, 5); 

    % run CKF to initialize measurements
    state0    = state_true(:, 1) + delta_state0;
    Nt_max    = round(3600 * fs);  % [sec]

    [P0_new, state0_new] = initialize_filter_CKF(time, state0, Y_true,...
         signal_err, R0, P0, Nt_max, planetParams, Cnm_filt, Snm_filt,...
         instrument_alig, NB_EARTH_mat, NB_MOON_mat, Q0, Qb, mask);

% %     % test initial error and uncertainty
% %     err0      = abs(state0_new - state_true(:, 1));
% %     sigma3    = 3.*sqrt(diag(P0_new));

    % run EKF
    X0 = state0_new; P0 = P0_new;

    [X_EKF, P_EKF] = EKF_process(time, planetParams, Y_true, signal_err,...
        X0, P0, instrument_alig, NB_EARTH_mat, NB_MOON_mat, Cnm_filt,...
        Snm_filt, Q0, Qb, R0, mask);
end

