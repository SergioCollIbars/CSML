clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Based on derivation from Rummel and Gelderen (1992), test
% the orthogonalityo fo the SH tensor and compute the sensivity matrix to
% posiiton errors.

syms GM r R positive
syms lambda phi

% coefficients vector
n_max  = 2; Ncs = 0;
for j = 2:n_max
    Ncs =  Ncs + 2*j + 1;
end

% SH tensor components
Tzz = 0; Tzx = 0; Tzy = 0;
Tyy = 0; Txx = 0; Txy = 0;
for j = 2:n_max
    n = j;
    for k = 0:j
        m = k;
        Pnm   = NormFactor(n, m) * assocLegendre(n, m);
        dPlm  = NormFactor(n, m) * Diff_assocLegendre(n,m);
        ddPlm = NormFactor(n, m) * Diff2_assocLegendre(n,m);

        valC = sym(['C_' num2str(n) '_' num2str(m)]);
        valS = sym(['S_' num2str(n) '_' num2str(m)]);

        Tzz = Tzz + GM / (r^3) * (R^n) / (r^n) * Pnm * ...  % radial - radial
            (cos(k * lambda) * valC + sin(k * lambda) * valS ) * (n+1) * (n+2);
        Tzx = Tzx - GM / (r^3) * (R^n) / (r^n) * dPlm * ...  % radial - normal
            (cos(k * lambda) * valC + sin(k * lambda) * valS ) * (n+2);
        Tzy = Tzy - GM / (r^3) * (R^n) / (r^n) * Pnm * ...  % radial - tangential
            (-sin(k * lambda) * valC + cos(k * lambda) * valS ) * (n+2) * m / cos(phi);
        Tyy = Tyy + GM / (r^3) * (R^n) / (r^n) * ...        % normal - normal
            (cos(k * lambda) * valC + sin(k * lambda) * valS ) * (ddPlm -(n+1)*Pnm);
        Txy = Txy + GM / (r^3) * (R^n) / (r^n) * ...        % normal - tangential
            (-sin(k * lambda) * valC + cos(k * lambda) * valS ) * m/cos(phi) * (ddPlm +cos(phi)*Pnm); 
    end
end

% form SH groups
Z1 = Tzz;
Z2 = [Tzx;Tzy];
Z3 = [Txx - Tyy; 2*Txy];
Z4 = [Z1; Z2];

% compute Gravity field partials
Z = Z1;
Nm = length(Z);
syms Hc [Nm Ncs]
syms Hr [Nm, 1]
count = 1;
for j = 2:n_max
    n = j;
    for k = 0:j
        m = k;

        valC = sym(['C_' num2str(n) '_' num2str(m)]);
        for m = 1:Nm
            Hc(m, count) = diff(Z(m, 1), valC);  % Partials w.r.t Cnm
        end
        count = count + 1;
    end
end
for j = 2:n_max
    n = j;
    for k = 1:j
        m = k;

        valS = sym(['S_' num2str(n) '_' num2str(m)]);

        for m = 1:Nm
            Hc(m, count) = diff(Z(m, 1), valS);  % Partials w.r.t Snm
        end
        count = count + 1;
    end
end
syms B;
Tprr = - diff(Tzz, "phi") + 2*r * (Tzx);
Trrr = [(-(3+n)/r) * Z1; (-(4+n)/r) * Z2];
Tlrr = Tzy * (cos(phi)*B - 2/r);
Hr(:, 1) = Tlrr;

% Compute visibility Matrix
Ax = Hc' * Hc; syms AxInt [Ncs, Ncs]
for row = 1:Ncs
    for col = 1:Ncs
        expression = Ax(row, col) * cos(phi);
        I1 = int(expression, phi, -pi/2, pi/2);
        I  = int(I1, lambda, 0, 2*pi);
        AxInt(row, col) = I;
    end
end

% Compute Mxc Matrix
Mxc = Hc' * Hr; syms MxInt [Ncs, 1]
for row = 1:Ncs
    for col = 1:1
        expression = Mxc(row, col) * cos(phi);
        I1 = int(expression, phi, -pi/2, pi/2);
        I  = int(I1, lambda, 0, 2*pi);
        MxInt(row, col) = I;
    end
end

% Sensitivityy matrix
Px = diag(diag(1./AxInt));
Sxc = Px * MxInt;

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