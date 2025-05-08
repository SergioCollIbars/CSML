clear;
clc;
close all;

syms GM r R 
syms C00 C20 C30 C40 C50
syms lambda phi

% coefficients vector
Cvec = [C00, C20, C30, C40, C50];

% compute Hc visibility matrix
syms Hc [1, 4]

% n = 0
Pnm = assocLegendre(0, 0);
val = 2 * GM /r^3 * Pnm;
Hc(1, 1) = val;

% n = 2 to n = 5
n = 2; m = 1;
for j = 2:5
    Pnm = assocLegendre(n, m);
    Pnm = NormFactor(n, m) * Pnm;
    dPlm = Diff_assocLegendre(n,m);
    dPlm =  NormFactor(n, m) * dPlm;
    val = (2+n)*(1+n)*GM*R^n/(r^(3+n)) * Pnm;
    Hc(1, j) = val;
    n = n + 1;
end

% compute Hr visibility matrix
syms Hr [1, 1]

% n = 0
Pnm = assocLegendre(0, 0);
val = -6 * GM /r^4 * Pnm * C00;
Hr(1, 1) = val;

% n =2 to n =5
n = 2;
for j = 2:5
    Pnm = assocLegendre(n, 1);
    dPlm = Diff_assocLegendre(l,m);
    %Pnm = NormFactor(n, 0) * Pnm;
    val = -(3+n)*(2+n)*(1+n)*GM*R^n/(r^(4+n)) * Pnm * Cvec(j);
    Hr(1, 1) = Hr(1, 1) + val;
    n = n +1;
end

% compute Ax matrix
syms AxInt [4, 4]
Ax = transpose(Hc) * Hc;
% intergrate over the whole planet surface
for j = 1:5
    val = Ax(j ,j) * cos(phi);
    expression = int(val, phi, -pi/2, pi/2);
    solution = int(expression, lambda, 0, 2*pi);
    AxInt(j, j) = solution;
end

Px = diag(1./diag(AxInt));

% compute Mxc matrix
Mxc  = transpose(Hc) * Hr;
syms MxcInt [4, 1]
% intergrate over the whole planet surface
for j = 1:5
    val = Mxc(j ,1)*cos(phi);
    expression = int(val, phi, -pi/2, pi/2);
    solution = int(expression, lambda, 0, 2*pi);
    MxcInt(j, 1) = solution;
end

% Correlation funciton
Sxc = -Px * MxcInt;
Pxc = Sxc * transpose(Sxc);
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