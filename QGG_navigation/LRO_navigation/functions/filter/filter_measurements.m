function [X_EKF, P_EKF, posfit,I_ALIG] = filter_measurements(metaData_file,...
    time, state_true, NB_EARTH_mat, NB_MOON_mat, Cnm_true, Snm_true, ...
    Y_true, signal_err, planetParams, instrumentParams_GG, instrumentParams_ST)
    %% FILTER GRADIOMETER MEASUREMENTS
    
    % load filter parameters
    [R0, P0, Q0, Qb, delta_state0, Cnm_filt, Snm_filt, I_ALIG, ...
        sigma_att, n_max] = ...
        loadFilterParams(metaData_file, instrumentParams_GG,...
        Cnm_true, Snm_true);
    planetParams_filter    = planetParams;
    planetParams_filter(5) =  n_max;
    
    % get measurement mask
    mask = instrumentParams_GG(:, 1);

    % create attitude error
    noise_att  = normrnd(0, instrumentParams_ST(1), ...
        [3, length(state_true(1, :))]);
    offset_att = instrumentParams_ST(2).*ones(3, length(time));

    % true SF
    SF        = state_true(13:18, :);

    % nominal state at initial time
    state0    = state_true(:, 1) + delta_state0;

    P0_new     = P0;
    state0_new = state0;

    % run EKF
    % % X0 = state0_new; P0 = P0_new;
    X0 = state0_new(1:12); P0 = P0_new(1:12, 1:12);

    [X_EKF, P_EKF, posfit] = EKF_process(time, planetParams_filter, ...
        Y_true, signal_err, SF, noise_att, offset_att, sigma_att, X0, P0, ...
        I_ALIG, NB_EARTH_mat, NB_MOON_mat, Cnm_filt, Snm_filt, Q0, Qb,...
        R0, mask);
end

