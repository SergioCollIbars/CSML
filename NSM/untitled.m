clear;
clc;

syms x phi lambda

%specify initial degree and final 
n_init = 8;
n_final = 30;
m = 6;

% initial
[~, Pnm] = assocLegendre(n_init,m, 1);
for n = m:n_final
    [~, Plm] = assocLegendre(n, m, 1);
    
    expr = Pnm * Plm / sqrt(1-(x^2));
% %     expr = Pnm * Plm;

    I = int(expr, 'x', -1, 1);
    disp(I)
end

n1 = 5;
m1 = 2;
[~, Pnm] = assocLegendre(n1,m1, 1);
% % [~, Pnm1] = assocLegendre(n1+1,m1+1, 0);
% % [~, Pnm2] = assocLegendre(n1+1,m1-1, 0);
% % 
% % Pnm3 = - 1/(2*m1) * (Pnm1 + (n1-m1+1)*(n1-m1+2)*Pnm2);

n2 = 5;
m2 = m1 + 1;

[~, Plk] = assocLegendre(n2,m2, 1);

expr = Pnm * Plk;
% % expr = Pnm3 * Plk;
% % expr = Pnm * Plk * cos(phi);

I = int(expr, 'x', -1, 1);
disp(I)


%% FUNCTIONS
function [plm_phi, plm_x] = assocLegendre(l,m, norm)
    % get symbolic associated legendre function P_lm(x) based on
    % legendre function P_l(x)
    
    syms x phi;
    
    % get symbolic form of Legendre function P_l(x)
    leg = legendreP(l, x);

    % differentiate it m times
    legDiff = diff(leg,x,m);

    % calculate associated legendre function P_lm(x)
    plm =  ((1 - x^2)^(m/2))*legDiff;
    
    plm_x = plm;
    plm_phi = subs(plm, x, sin(phi));
    
    if(norm)
        a = (l + 1/2)*factorial(l-m);
        b = factorial(l+m);
    
    % %     a = 2*factorial(l-m);
    % %     b = factorial(l + m);
        N = (-1)^m * sqrt(a/b);
        plm_x = N * plm_x;
    end
end