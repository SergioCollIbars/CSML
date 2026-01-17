function [Y, H, Hr] = compute_measurements_filter(planetParams, ...
                          time, state, Cnm, Snm, BN_mat, BN_MOON_mat)
    % extract params
    Nt = length(time);

    GM_M        = planetParams(2); 
    R_M         = planetParams(4);
    n_max       = planetParams(5);
    normalized  = planetParams(6);

    %% compute measurements
    k = 1;  % time index
    Y = ones(6, Nt) * NaN;
    
    maxInd = 3 * k; minInd = maxInd - 2;
    BODYMOON_J2000 = BN_MOON_mat(minInd:maxInd, :);

    r = state(k, 1:3)'; v = state(k, 4:6)';
    [Y_J2000, ~] = gradiometer_meas(time(k) ,[GM_M, R_M, n_max, normalized],...
        BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm{2}, Snm{2});

    % rotation to Instrument frame
    BN = BN_mat(minInd:maxInd, :);

    T_ACI = [Y_J2000(1),Y_J2000(2),Y_J2000(3);...
             Y_J2000(4),Y_J2000(5),Y_J2000(6);...
             Y_J2000(7),Y_J2000(8),Y_J2000(9)];
    T_B   = BN * T_ACI * BN';

    Y(:, k) = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);...
               T_B(2,3);T_B(3,3)]./1E-12;   % [mE]

    %% compute measurement partials
    BODYMOON_BODY = BODYMOON_J2000 * BN'; r_ACI = r;
    [Hpos] = compute_posPartials(n_max, normalized, Cnm{2}, Snm{2}, R_M, GM_M, ...
                r_ACI, BODYMOON_J2000, BODYMOON_BODY);
    Hp     = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;
    H      = [Hp, zeros(6, 3), eye(6)];


    %% compute attitude partials
    [Hrot] = compute_rotPartials_analy(Y_J2000, BN);
    Hr     = [Hrot(1:3, :); Hrot(5:6, :);Hrot(9, :)]./1E-12;
end

