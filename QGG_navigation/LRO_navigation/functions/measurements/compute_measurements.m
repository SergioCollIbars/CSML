function [Y, bias] = compute_measurements(instrumentParams, planetParams, ...
    time, state, Cnm, Snm, BN_mat, NB_EARTH_mat, NB_MOON_mat)
    
    % extract params
    Nt = length(time);

    GM_M        = planetParams(2); 
    R_M         = planetParams(4);
    n_max       = planetParams(5);
    normalized  = planetParams(6);

    % generate white noise
    sigma_m = sqrt(instrumentParams(1, 5)) * instrumentParams(1, 2);
    noise_white = normrnd(0, sigma_m, [6, Nt]);

    % generate flicker (alpha = 2) noise
    S_w        = instrumentParams(1, 2);                           % [mE/sqrt(Hz)]

    fs         =  instrumentParams(1, 5);                          % [Hz]
    f_min      =  instrumentParams(1, 6);                          % [Hz]

    sigma_eta2 = 4 * S_w^2 * sin(pi * f_min / fs)^2;

    eta             = sqrt(sigma_eta2) * randn(6, Nt);    % RW increments
    noise_flicker   = cumsum(eta')';                      % random walk

    %% compute measurements
    Y = ones(6, Nt) * NaN;
    for k = 1:Nt
        maxInd = 3 * k; minInd = maxInd - 2;
        BODYMOON_J2000 = NB_MOON_mat(minInd:maxInd, :)';
    
        r = state(k, 1:3)'; v = state(k, 4:6)';                             % WARNING: Only accounting for Moon gravity gradient
        [Y_J2000, ~] = gradiometer_meas(time(k) ,...
            [GM_M, R_M, n_max, normalized],...
            BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm{2}, Snm{2});

        % rotation to Instrument frame
        maxInd = 3 * k; minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);

        T_ACI = [Y_J2000(1),Y_J2000(2),Y_J2000(3);...
                 Y_J2000(4),Y_J2000(5),Y_J2000(6);...
                 Y_J2000(7),Y_J2000(8),Y_J2000(9)];
        T_B   = BN * T_ACI * BN';

        Y(:, k) = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);...
                   T_B(2,3);T_B(3,3)]./1E-12;   % [mE]
    end
    
    % add noise
    cnst_bias = instrumentParams(:, 3).*ones(6, Nt);
    Y         = Y + noise_white + cnst_bias + noise_flicker;
    
    % update true bias
    bias = noise_flicker' + cnst_bias';
end