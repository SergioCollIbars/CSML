%% One-shot J2 disturbing potential per unit mass:
% Kaula (cos Psi with u = omega + f) vs Classical (r,phi,lambda)

clear; clc;

% ---- Body constants (Earth defaults; change as needed) ----
mu  = 3.986004418e5;    % [km^3/s^2]
R   = 6378.137;         % [km]
J2  = 1.08262668e-3;    % [-]

% ---- Choose ONE set of orbital elements (radians) ----
a     = 7000;                     % [km]
e     = 0.1;                        % set 0 for circular test
i     = deg2rad(90);              % inclination [rad]
Omega = deg2rad(40);              % RAAN [rad]
omega = deg2rad(30);              % arg. perigee [rad]
M     = deg2rad(20);              % mean anomaly [rad]
eng   = -mu / (a);              % specific energy 

% ---- (Optional) body sidereal angle to go ECI -> body-fixed (for lambda)
theta = deg2rad(20);               % set to 0 if you like

% ---- Solve Kepler, get f, r, u ----
E  = keplerE(M, e);
f  = 2*atan2( sqrt(1+e).*sin(E/2), sqrt(1-e).*cos(E/2) );
r  = a*(1 - e*cos(E));
u  = omega + f;                   % argument of latitude

% ---- Kaula-style J2 via cos(Psi) with l=2, m=0, q=0 ----
da_dr = a / r;

[R_J2_Kaula] = computePotential_kaula(mu, R, J2, 2, 0, a, i, e, ...
    omega, Omega, M, theta);
dR_J2_dr_Kaula = (-(3) / a * R_J2_Kaula) * da_dr;
ddR_J2_ddr_Kaula = ((12) / (a^2) * R_J2_Kaula) * (da_dr)^2;


% ---- Classical J2 via (r, phi, lambda) and Legendre P2 ----
% Build ECI position from elements, rotate to body-fixed to get lambda
[rECI, ~] = coe2rv(a,e,i,Omega,omega,M,mu);
rmag = norm(rECI);
R3th = rotz(-theta);
rBF  = R3th * rECI;
phi  = asin(rBF(3)/rmag);         % geocentric latitude
lambda = atan2(rBF(2), rBF(1));

R_J2_classical = (mu*J2*R^2)/(2*rmag^3) * ( 3*sin(phi)^2 - 1 );
dRJ2_dr_class = -3*(mu*J2*R^2)/(2*rmag^4) * ( 3*(sin(phi)^2) - 1 );
ddRJ2_ddr_class = 12*(mu*J2*R^2)/(2*rmag^5) * ( 3*(sin(phi)^2) - 1 );

% ---- Compare ----
abs_err = abs(R_J2_Kaula - R_J2_classical);
rel_err = abs_err / max(1e-16, abs(R_J2_classical));
fprintf('J2 potential check:\n  Kaula  : %.12e km^2/s^2\n  Classic: %.12e km^2/s^2\n', ...
        R_J2_Kaula, R_J2_classical);
fprintf('  |abs diff| = %.3e,   rel = %.3e\n', abs_err, rel_err);

%% ---------- helpers ----------
function E = keplerE(M,e)
% Solve Kepler's equation M = E - e sin E (scalar)
M = mod(M, 2*pi);
if e < 0.8, E = M; else, E = pi; end
for k=1:50
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f/fp;
    E  = E + dE;
    if abs(dE) < 1e-14, break; end
end
end

function [rI, vI] = coe2rv(a,e,i,Omega,omega,M,mu)
% Elements -> ECI r,v (scalar)
E = keplerE(M,e);
f = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
p = a*(1 - e^2);
r = p/(1 + e*cos(f));
r_pf = [r*cos(f); r*sin(f); 0];
v_pf = sqrt(mu/p)*[-sin(f); e+cos(f); 0];

Q = rotz(Omega)*rotx(i)*rotz(omega);
rI = Q*r_pf; vI = Q*v_pf;
end

function R = rotx(a)
ca=cos(a); sa=sin(a);
R = [1 0 0; 0 ca -sa; 0 sa ca];
end
function R = rotz(a)
ca=cos(a); sa=sin(a);
R = [ca -sa 0; sa ca 0; 0 0 1];
end

function [U] = computePotential_kaula(GM, Ref, J2, l, m, a, i, e, ...
    omega, Omega, M, theta)
    Amplitude = GM * Ref^l / ((a)^(l+1));
    sum_over_p = 0;
    for p = 0:l
        % inclination function
        F = ffun2(l, m, p, i);

        % eccentircity function
        sum_over_q = 0;
        for q = -20:20
            G = gfun2(l,p,q,e,10);
            Psi = (l - 2*p)*omega + (l - 2*p + q)*M + m*(Omega - theta) ;
            sum_over_q = sum_over_q + J2 * cos(Psi) * G;
        end
        sum_over_p = sum_over_p + sum_over_q * F;
    end

    U = Amplitude * sum_over_p;
end


function f_num = ffun2(l, m, p, i)
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

function C = nchoosek_safe(n, r)
% binomial with automatic 0 outside range, using log-factorials for stability
    if r < 0 || n < 0 || r > n
        C = 0.0;
    else
        C = exp(gammaln(n + 1) - gammaln(r + 1) - gammaln(n - r + 1));
    end
end

function g_num=gfun(l,p,q, ecc)

syms e k rq rp d s t

if p>l/2
    pp=l-p;
    qp=-q;
else
    pp=p;
    qp=q;
end


if qp<0
    hp=k;
    hq=k-qp;
else
    hp=k+qp;
    hq=k;
end



beta=e/(1+sqrt(1-e^2));

g=(-1)^abs(q)*(1+beta^2)^l*beta^abs(q)*symsum(symsum(nchoosek(2*pp-2*l,hp-rp)*(((-1)^rp)/factorial(rp))*((l-2*pp+qp)*e/(2*beta))^rp,rp,0,hp)*symsum(nchoosek(-2*pp,hq-rq)/factorial(rq)*((l-2*pp+qp)*e/(2*beta))^rq,rq,0,hq)*beta^(2*k),k,0,10);

g_num = double(subs(g, 'e', ecc));
end

function g_num = gfun2(l,p,q,ecc,Kmax)
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

% ---------- Helper: generalized binomial for integer n, k >= 0 ----------
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
