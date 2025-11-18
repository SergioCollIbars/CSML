function [g_num] = compute_eccentricityFunction(l, p, q, ecc, Kmax)
%GFUN Numerical evaluation of the truncated series for g(l,p,q,e)
%   g_num = gfun(l,p,q,ecc)        % uses Kmax = 10 (to match your code)
%   g_num = gfun(l,p,q,ecc,Kmax)   % custom truncation on k-sum
%
% Inputs:
%   l,p,q : integers (as in your original)
%   ecc   : eccentricity (0 <= ecc < 1 recommended)
%   Kmax  : nonnegative integer (default 10)
%
% Notes:
% - Implements the same logic as your symbolic code:
%       pp, qp selection
%       hp,hq depend on k and sign(qp)
%       beta = e/(1+sqrt(1-e^2))
%       Double inner sums over rp and rq, outer sum over k=0..Kmax
% - Uses generalized binomial for integer n (possibly negative).
% - Factorials via gammaln for stability.

    if nargin < 5
        Kmax = 10;
    end

    % ---- pp, qp branches (same as your code) ----
    if p > l/2
        pp = l - p;
        qp = -q;
    else
        pp = p;
        qp = q;
    end

    % ---- beta and ratio (handle ecc=0 safely) ----
    if ecc == 0
        beta  = 0;
        ratio = 1;  % limit of e/(2*beta) as e->0 is 1
    else
        beta  = ecc/(1 + sqrt(1 - ecc^2));
        ratio = ecc/(2*beta);
    end

    pref_sign = (-1)^(abs(q));
    pref_mag  = (1 + beta^2)^l * beta^(abs(q));
    pref      = pref_sign * pref_mag;

    % alpha = ((l - 2*pp + qp) * e / (2*beta)) with safe ratio
    alpha = (l - 2*pp + qp) * ratio;

    % ---- Outer sum over k ----
    g_sum = 0.0;
    for k = 0:Kmax
        % hp, hq depend on sign of qp (exactly like your code)
        if qp < 0
            hp = k;
            hq = k - qp;  % note: qp<0 => hq >= k
        else
            hp = k + qp;
            hq = k;
        end

        % ---- Inner sums over rp=0..hp and rq=0..hq ----
        S1 = 0.0;
        for rp = 0:hp
            % binom(2*pp - 2*l, hp - rp) * (-1)^rp / rp! * alpha^rp
            b1   = binom_int(2*pp - 2*l, hp - rp);
            term = b1 * ((-1)^rp) * exp(-gammaln(rp + 1)) * (alpha^rp);
            S1   = S1 + term;
        end

        S2 = 0.0;
        for rq = 0:hq
            % binom(-2*pp, hq - rq) * 1/rq! * alpha^rq
            b2   = binom_int(-2*pp, hq - rq);
            term = b2 * exp(-gammaln(rq + 1)) * (alpha^rq);
            S2   = S2 + term;
        end

        % beta^(2k) (handle beta=0 cleanly: 0^(0)=1 for k=0; 0 for k>0)
        if beta == 0
            beta_pow = double(k == 0);
        else
            beta_pow = beta^(2*k);
        end

        g_sum = g_sum + (S1 * S2 * beta_pow);
    end

    g_num = pref * g_sum;
end

%% Auxiliary functions
function c = binom_int(n,k)
%BINOM_INT Generalized binomial coefficient for integer n (can be negative) and k>=0.
% Uses identity for negative n: C(n,k) = (-1)^k * C(k - n - 1, k)
    if k < 0 || floor(k) ~= k
        c = 0; return;
    end
    n = round(n); % ensure integer
    if n >= 0
        if k > n
            c = 0;
        else
            % standard integer binomial via gamma (exact for integers)
            c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));
        end
    else
        % negative n case
        c = ((-1)^k) * exp(gammaln(k - n) - gammaln(k+1) - gammaln(-n));
    end
end

