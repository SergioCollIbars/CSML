function [f_num] = compute_inclinationFunction(l, m, p, i)
% FFUN  Numeric Kaula inclination function F_{lmp}(i) using your sum formula.
% Input:  l,m,p integers with 0<=m<=l, 0<=p<=l; i in radians
% Output: f_num = F_{lmp}(i)

    % basic input checks
    if l < 0 || m < 0 || p < 0 || m > l || p > l
        error('ffun: bad indices. Need 0<=m<=l, 0<=p<=l.');
    end

    k   = floor((l - m)/2);   % your k = fix((l-m)/2)
    si  = sin(i);
    ci  = cos(i);
    f_num = 0.0;

    % outer sum over t = 0 .. min(p,k)
    tmax = min(p, k);
    for t = 0:tmax
        % guard: l - m - 2t must be >= 0 for factorials/powers to make sense
        lm2t = l - m - 2*t;
        if lm2t < 0
            continue
        end

        % prefactor:
        % ((2l-2t)!)/( t! (l-t)! (l-m-2t)! 2^(2l-2t) ) * sin(i)^(l-m-2t)
        logA = gammaln(2*l - 2*t + 1) ...
             - (gammaln(t + 1) + gammaln(l - t + 1) + gammaln(lm2t + 1)) ...
             - (2*l - 2*t) * log(2);
        A = exp(logA) * (si^lm2t);

        % middle sum over s = 0..m
        sum_s = 0.0;
        for s = 0:m
            comb_ms = nchoosek_safe(m, s);
            if comb_ms == 0, continue; end

            % inner sum over c
            % Your original had c=0..5, but the natural finite bounds are:
            %   0 <= c <= min(l-m-2t + s, p - t)
            cmax = min(lm2t + s, p - t);
            if cmax < 0, continue; end

            sum_c = 0.0;
            for c = 0:cmax
                term_c = nchoosek_safe(lm2t + s, c) * ...
                         nchoosek_safe(m - s, p - t - c) * ...
                         (-1)^(c - k);
                sum_c = sum_c + term_c;
            end

            sum_s = sum_s + comb_ms * (ci^s) * sum_c;
        end

        f_num = f_num + A * sum_s;
    end
end


%% Auxiliary functions
function C = nchoosek_safe(n, r)
% binomial with automatic 0 outside range, using log-factorials for stability
    if r < 0 || n < 0 || r > n
        C = 0.0;
    else
        C = exp(gammaln(n + 1) - gammaln(r + 1) - gammaln(n - r + 1));
    end
end
