function [Hpos] = compute_rotPartials(n_max, normalized, Cmat, Smat, Re, GM, r, ACAF_ACI)
    % output value
    Hpos = ones(9, 3) * NaN;
    ACI_ACAF = ACAF_ACI';

    eps = 1E-6;
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;

        Atpos = At./2;
        Atneg = - At./2; 

        [Rpos] = rotationMatrix(Atpos(1), Atpos(2), Atpos(3), [3, 2, 1]);
        [Rneg] = rotationMatrix(Atneg(1), Atneg(2), Atneg(3), [3, 2, 1]);

        [~, ~, ddUpos] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*r, Re, GM, ...
                                                normalized);
        ddUpos = ACI_ACAF * ddUpos * ACI_ACAF';
        ddUpos = Rpos' * ddUpos * Rpos;

        [~, ~, ddUneg] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*r, Re, GM, ...
                                                normalized);
        ddUneg = ACI_ACAF * ddUneg * ACI_ACAF';
        ddUneg = Rneg' * ddUneg * Rneg;

        H = (ddUpos - ddUneg)./(vecnorm(Atpos-Atneg));

        Hpos(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
            H(3,1); H(3,2);H(3,3)];
    end
end

