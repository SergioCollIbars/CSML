%% J2 disturbing potential per unit mass: Kaula in M (fixed ω via u−f) vs Classical
clear; clc;

% --- Body constants (Earth) ---
mu  = 3.986004418e5;     % km^3/s^2
R   = 6378.137;          % km
J2  = 1.08262668e-3;     % -

% --- Example ECI state ---
rECI = [6800; 100; 200];     % km
vECI = [ -1.0; 7.4; 1.0];    % km/s

% ---- Planet rotation angle (e.g., GMST/ERA in radians). Set as needed.
theta = deg2rad(37.0);  % put your GMST/ERA here

% ===== From (r,v): geometry, elements, and *consistent* angles =====
r = norm(rECI); khat=[0;0;1];
h = cross(rECI,vECI); hhat = h/norm(h);
n = cross(khat,h); nrm = norm(n);
if nrm<1e-14, e1=[1;0;0]; else, e1 = n/nrm; end        % toward ascending node
e2 = cross(hhat,e1);                                   % 90° ahead in plane

% Argument of latitude directly from r
u = atan2( dot(rECI,e2), dot(rECI,e1) );

% Classical elements (for a,e,f,M only)
[a,e,~,Omega,~,f,E,M] = rv2coe_basic(rECI, vECI, mu);

% Enforce ω so that u = ω + f holds numerically
omega = wrapTo2Pi(u - f);

% Inclination from h
i = acos(hhat(3));

rECEF = R3(theta) * rECI;                  % ECI -> ECEF (rotation about z)
rmag  = norm(rECEF);
phi   = asin(rECEF(3)/rmag);               % geocentric latitude from ECEF
% ===== Classical J2 =====
R_J2_class = (mu*J2*R^2)/(2*r^3) * ( 3*sin(phi)^2 - 1 );
dRJ2_dr_class = -3*(mu*J2*R^2)/(2*r^4) * ( 3*(sin(phi)^2) - 1 );

% ===== Kaula J2 =====
[R_J2_kaulaM] = computePotential_kaula(mu, R, J2, 2, 0, a, i, e, ...
    omega, Omega, M, 0);
da_dr = a / r;
dR_J2_dr_Kaula = (-(3) / a * R_J2_kaulaM) * da_dr;

% ===== Consistency diagnostics (should be ~machine-eps) =====
fprintf('check 1: z/r - sin(i)sin(u) = %.3e\n', rECI(3)/r - sin(i)*sin(u));
cu = cos(2*u);
cu_from_of = cos(2*omega)*cos(2*f) - sin(2*omega)*sin(2*f);
fprintf('check 2: cos2u(id) - cos2u(of) = %.3e\n', cu - cu_from_of);

% ===== Compare =====
abs_diff = abs(R_J2_kaulaM - R_J2_class);
rel_diff = abs_diff / max(1e-16, abs(R_J2_class));
fprintf('J2 per-unit-mass:\n  Kaula(M):   %.12e km^2/s^2\n  Classical:  %.12e km^2/s^2\n', ...
        R_J2_kaulaM, R_J2_class);
fprintf('  |abs diff|=%.3e   (rel=%.3e)\n', abs_diff, rel_diff);
abs_diff = abs(dR_J2_dr_Kaula - dRJ2_dr_class);
fprintf('J2 acceleration per-unit-mass:\n  Kaula(M):   %.12e km^2/s^2\n  Classical:  %.12e km^2/s^2\n', ...
        dR_J2_dr_Kaula, dRJ2_dr_class);
fprintf('  |abs diff|=%.3e   (rel=%.3e)\n', abs_diff, rel_diff);


% ---------- helpers ----------
function R = R3(a), ca=cos(a); sa=sin(a); R=[ca -sa 0; sa ca 0; 0 0 1]; end

function [a,e,i,RAAN,omega,f,E,M] = rv2coe_basic(rI,vI,mu)
r=norm(rI); v=norm(vI);
h=cross(rI,vI); h2=dot(h,h); hmag=sqrt(h2);
k=[0;0;1]; n=cross(k,h); nmag=norm(n);
evec=(1/mu)*((v^2-mu/r)*rI - dot(rI,vI)*vI); e=norm(evec);
eps=v^2/2 - mu/r; a=-mu/(2*eps);
i=acos(h(3)/hmag);
RAAN=atan2(n(2),n(1)); if RAAN<0, RAAN=RAAN+2*pi; end
omega=atan2( dot(cross(n,evec),h)/(hmag*nmag), dot(n,evec)/(nmag*e) );
if omega<0, omega=omega+2*pi; end
f=atan2( dot(cross(evec,rI),h)/(hmag*e*r), dot(evec,rI)/(e*r) ); if f<0, f=f+2*pi; end
E=2*atan2( sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2) ); if E<0, E=E+2*pi; end
M=E - e*sin(E); M=mod(M,2*pi);
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

