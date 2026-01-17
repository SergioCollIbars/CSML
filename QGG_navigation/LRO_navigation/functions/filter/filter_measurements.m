function [X_EKF, P_EKF, posfit,I_ALIG] = filter_measurements(metaData_file,...
    time, state_true, NB_EARTH_mat, NB_MOON_mat, Cnm_true, Snm_true, ...
    Y_true, signal_err, planetParams, instrumentParams_GG, instrumentParams_ST)
    %% FILTER GRADIOMETER MEASUREMENTS
    
    % load filter parameters
    [R0, P0, Q0, Qb, delta_state0, Cnm_filt, Snm_filt, I_ALIG, sigma_att] = ...
        loadFilterParams(metaData_file, planetParams, instrumentParams_GG,...
        Cnm_true, Snm_true);
    
    % get measurement mask
    mask = instrumentParams_GG(:, 1);
% %     fs   = instrumentParams_GG(1, 5); 

    % create attitude error
    noise_att  = normrnd(0, instrumentParams_ST(1), ...
        [3, length(state_true(1, :))]);
    offset_att = instrumentParams_ST(2).*ones(3, length(time));

    % nominal state at initial time
    state0    = state_true(:, 1) + delta_state0;

% %     % run CKF to initialize measurements
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
        noise_att, offset_att, sigma_att, X0, P0, I_ALIG, NB_EARTH_mat, ...
        NB_MOON_mat, Cnm_filt, Snm_filt, Q0, Qb, R0, mask);
end

