clear;
clc;
close all;

syms phi;

% specify max k 
k_max  = 150;

% specify n range
n_range = 0:8;

sol = ones(length(n_range), length(n_range))*NaN;
for j = 1:length(n_range)
    n = n_range(j);
    disp('Computing degree .... ' + string(n))
    for m = 0:n
        % compute legendre functions
        Pnm   = NormFactor(n, m) * assocLegendre(n, m);
        
        % compute expression
        if(m == 0)
            F = 2 *  NormFactor(0, m) * Diff_assocLegendre(0,m);
            for k  = 1:k_max
                dPlm  = NormFactor(k, m) * Diff_assocLegendre(k,m);
                F = F  + dPlm * (k+2) * ((2*k + 1)^(0.5))/ (k^2);
            end
        else
            F = 0;
            for k  = m:k_max
                dPlm  = NormFactor(k, m) * Diff_assocLegendre(k,m);
                F = F  + dPlm * (k+2) * ((2*k + 1)^(0.5))/ (k^2);
            end
        end
    
        % compute final expression
        I = Pnm * F * cos(phi);
        
        sol(m+1, j) = int(I, phi, -pi/2, pi/2);
    end
end

figure()
plot(n_range, sol, 'LineWidth', 2)

%%  FUNCTIONS
function plm = assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
end

function dPlm = Diff_assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
    dPlm = diff(plm, phi);
end

function ddPlm = Diff2_assocLegendre(l,m)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l,x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm = ((1 - x^2)^(m/2))*legDiff;

    plm = subs(plm, x, sin(phi));
    dPlm = diff(plm, phi);
    ddPlm = diff(dPlm, phi);
end

function [N] = NormFactor(n, m)
    % Description: given degree, n and order, m compute the normalice
    % factor
    if(m == 0)
        delta = 1;
    else
        delta = 0;
    end
    fac1 = factorial(n - m);
    fac2 = factorial(n + m);
    N = ((2 - delta)*(2*n + 1) * fac1 /fac2)^(0.5);
end