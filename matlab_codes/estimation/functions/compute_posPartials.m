function [Hpos] = compute_posPartials(n_max, normalized, Cmat, Smat, Re, GM, r, ACAF_ACI, ACAF_B)
    % output value
    Hpos = ones(9, 3) * NaN;
    ACI_ACAF = ACAF_ACI';
    B_ACI = ACAF_B' * ACAF_ACI;

    eps = 1E-6;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = r + Ar./2;       % [ACI]
        rneg = r - Ar./2;       % [ACI]

        rposB = B_ACI * rpos;   % [Body]
        rnegB = B_ACI * rneg;   % [Body]

        [~, ~, ddUpos] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*rpos, Re, GM, ...
                                                normalized);
        ddUpos = ACAF_B' * ddUpos * ACAF_B;
        [~, ~, ddUneg] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*rneg, Re, GM, ...
                                                normalized);
        ddUneg = ACAF_B' * ddUneg * ACAF_B;
        H = (ddUpos - ddUneg)./(vecnorm(rposB-rnegB)); % position partials in body frame

        Hpos(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
            H(3,1); H(3,2);H(3,3)];
    end
end

