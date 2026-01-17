function [X_EKF, P_EKF, posfit,I_ALIG] = filter_measurements(metaData_file,...
    time, state_true, NB_EARTH_mat, NB_MOON_mat, Cnm_true, Snm_true, ...
    Y_true, signal_err, planetParams, instrumentParams)
    %% FILTER GRADIOMETER MEASUREMENTS
    
    % load filter parameters
    [R0, P0, Q0, Qb, delta_state0, Cnm_filt, Snm_filt, I_ALIG, err_att] = ...
        loadFilterParams(metaData_file, planetParams, instrumentParams,...
        Cnm_true, Snm_true);
    
    % get measurement mask
    mask = instrumentParams(:, 1);
    fs   = instrumentParams(1, 5); 

    % create attitude error
    sigma_att  = err_att(1); bias_att = err_att(2);
    noise_att  = normrnd(0, sigma_att, [3, length(state_true(1, :))]);
    offset_att = bias_att.*ones(3, length(time));

% %     noise_att = normrnd(0, sigma_att, [3, length(state_true(1, :))]);

% %     % run CKF to initialize measurements
    state0    = state_true(:, 1) + delta_state0;
% %     Nt_max    = round(3600 * fs);  % [sec]
% % 
% %     [P0_new, state0_new] = initialize_filter_CKF(time, state0, Y_true,...
% %          signal_err, noise_att, R0, P0, Nt_max, planetParams, Cnm_filt, Snm_filt,...
% %          I_ALIG, NB_EARTH_mat, NB_MOON_mat, Q0, Qb, mask);
% % 
% %     % test initial error and uncertainty
% %     err0      = abs(state0_new - state_true(:, 1));
% %     sigma3    = 3.*sqrt(diag(P0_new));

    P0_new     = P0;
    state0_new = state0;

    % run EKF
    X0 = state0_new; P0 = P0_new;

    [X_EKF, P_EKF, posfit] = EKF_process(time, planetParams, Y_true, signal_err,...
        noise_att, offset_att, sigma_att, X0, P0, I_ALIG, NB_EARTH_mat, NB_MOON_mat, Cnm_filt,...
        Snm_filt, Q0, Qb, R0, mask);
end

