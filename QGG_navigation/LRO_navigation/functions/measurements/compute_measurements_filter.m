function [Y, H, Hc, Hp] = compute_measurements_filter(planetParams, ...
    poleParams, time, state, Cnm, Snm, BN_mat)
    % extract params
    Nt = length(time);

    GM          = planetParams(1); 
    Re          = planetParams(2);
    n_max       = planetParams(3);
    normalized  = planetParams(4);

    W           = poleParams(1);
    W0          = poleParams(2);
    RA          = poleParams(3);
    DEC         = poleParams(4);

    %% compute measurements
    k = 1;  % time index
    Y = ones(6, Nt) * NaN;

    Wt = W0 + W * time(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    r = state(k, 1:3)'; v = state(k, 4:6)';
    [Y_ACI, ~] = gradiometer_meas(time(k) ,[GM, Re, n_max, normalized],...
        ACAF_ACI, [r', v'], zeros(9, 1), Cnm, Snm);

    % rotation to Instrument frame
    maxInd = 3 * k; minInd = maxInd - 2;
    BN = BN_mat(minInd:maxInd, :);

    T_ACI = reshape(Y_ACI, [3, 3]);
    T_B   = BN * T_ACI * BN';

    Y(:, k) = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);T_B(2,3);T_B(3,3)]./1E-12 + ...
        state(k, 7:12)';

    %% compute measurement partials
    ACAF_BODY = ACAF_ACI * BN'; r_ACI = r;
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                r_ACI, ACAF_ACI, ACAF_BODY);
    Hp     = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;
    H      = [Hp, zeros(6, 3), eye(6)];

    %% Compute measurement partials
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                r_ACI, ACAF_ACI, ACAF_BODY);
     Hp   = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;
     
    %% compute consider covariance measurement partials
    r_ACAF   = ACAF_ACI * r_ACI;
    [~, Hnm] = potentialGradient_Cnm(n_max, r_ACAF, Re, GM, ...
                        ACAF_BODY', normalized);
    Hc       = [Hnm(1:3, :); Hnm(5:6, :); Hnm(9, :)]./1E-12;
end

