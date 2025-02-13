clear;
clc;
close all
%%          LEGENDRE POLINOMIALS TEST

syms lambda phi n
n1 = 2;
m1 = 0;

n2 = 2;
m2 = 0;

normalized = 1;

P1 = assocLegendre(n1, m1);
P2 = assocLegendre(n2, m2);

dP1 = Diff_assocLegendre(n1, m1);
dP2 = Diff_assocLegendre(n2, m2);

ddP1 = Diff2_assocLegendre(n1, m1);
ddP2 = Diff2_assocLegendre(n2, m2);

if(normalized)
    P1 = NormFactor(n1, m1) * P1;
    P2 = NormFactor(n2, m2) * P2;
    dP1 = NormFactor(n1, m1) * dP1;
    dP2 = NormFactor(n2, m2) * dP2; 
    ddP1 = NormFactor(n1, m1) * ddP1;
    ddP2 = NormFactor(n2, m1) * ddP2;
end

% create hi functions (meas functions)
h1 = (n + 1) * (n + 2) * P1 * cos(m1 * lambda);
h2 = -dP1 * (n + 2) * cos(m1 * lambda);
h4 = (ddP1 - (n + 1) * P1) * cos(m1 * lambda);
h6 = (dP1/tan(phi) - P1*(n+1)) * cos(m1 * lambda) * 0;
H = (h1)^2 + (h2)^2 + (h4)^2 + (h6)^2;
H = H * cos(phi);
expression = int(H, phi, -pi/2, pi/2);
solution = int(expression, lambda, 0, 2*pi);
disp('Solution H: ' + string(simplify(solution)));

% compare subs solution and polynomial
subSol = subs(solution, n, n1);
pol = (n1+1)^2 *(n1+2)^2 + (n1+2)^2*(n1+1)*n1 + ...
    (n1+1)*n1*(n1^2 -0.5) + 2*n1^2*(n1+1) + (n1+1)^2;
pol = pol*4*pi;

% Compute Geodesic functions
Y1 = P1 * cos(m1 * lambda);
Y2 = P2 * cos(m2 * lambda);

dY1 = dP1 * cos(m1 * lambda);
dY2 = dP2 * cos(m2 * lambda);

ddY1 = ddP1 * cos(m1 * lambda);
ddY2 = ddP2 * cos(m2 * lambda);

expression = Y1 * Y2 * cos(phi);

expression = int(expression, phi, -pi/2, pi/2);
solution = int(expression, lambda, 0, 2*pi);
disp('Solution Y: ' + string(solution));

expression = dY1 * dY2 * cos(phi);

expression = int(expression, phi, -pi/2, pi/2);
solution = int(expression, lambda, 0, 2*pi);
disp('solution dY: ' + string(solution));

expression = ddY1 * ddY2 * cos(phi);

expression = int(expression, phi, -pi/2, pi/2);
solution = int(expression, lambda, 0, 2*pi);
disp('solution ddY: ' + string(solution));

expression = (ddY1 * Y2) * cos(phi);
expression = int(expression, phi, -pi/2, pi/2);
solution = int(expression, lambda, 0, 2*pi);
disp('solution special expression: ' + string(solution));



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