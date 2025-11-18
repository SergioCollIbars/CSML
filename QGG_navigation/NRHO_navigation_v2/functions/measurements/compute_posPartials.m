function [H] = compute_posPartials(planetParams, BN, C_mat, S_mat, t, x, posE, posM, posS)
    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [ACI]
        rneg = x - Ar./2;   % [ACI]
        [ddU_pos_N] = compute_nBody(rpos ,t, C_mat, S_mat, planetParams, posE, posM, posS);
        [ddU_neg_N] = compute_nBody(rneg ,t, C_mat, S_mat, planetParams, posE, posM, posS);
        
        ddU_pos = BN * ddU_pos_N * BN';
        ddU_neg = BN * ddU_neg_N * BN';
        Ht = (ddU_pos - ddU_neg)./(vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end

