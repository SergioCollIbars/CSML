function [R] = compute_whitening(N, sigma_W, sigma_RW)
    R = zeros(N,N);
    qb = sigma_RW^2;
    for k = 1:N
        for l = 1:N
            R(k,l) = (sigma_W^2)*(k==l) + qb*min(k,l);
        end
    end
end

