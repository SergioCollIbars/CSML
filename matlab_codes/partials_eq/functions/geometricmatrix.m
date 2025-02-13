function [GEOM] = geometricmatrix(n_max, r, phic, lambda, GM, R, normalized)
    % Description: computes the geometric matrix for a given set of SH
    % harmonics degree & order, map coordinates and planet parameters.

    % # of states
    Nc = - 2;                              
    for j = 1:n_max+1
        Nc = Nc + j;
    end
    Ns = Nc - n_max;
    Nx = Nc + Ns;
    
    % define Geometric matrix
    GEOM = zeros(6, Nx);
    
    % point mass value. n = m = 0
    n = 0;
    m = 0;
    c = cos(m * lambda);
    s  = sin(m * lambda); 
    sP = sin(phic);
    cP = cos(phic);
    tP = sP/cP;
    ctP = 1/tP;
    P = double(subs(assocLegendre(n, m), phic));
    dP = double(subs(Diff_assocLegendre(n, m), phic));
    ddP = double(subs(Diff2_assocLegendre(n, m), phic));

    if(normalized)
        P = NormFactor(n, m) * P;
        dP = NormFactor(n, m) * dP;
        ddP = NormFactor(n, m) * ddP;
    end

    Bn = (GM/r^3) * (R/r)^n;

    GEOM(:, 1) = ...
    [(n + 1)*(n + 2)*P*c;...
    -dP*(n+2)*c;...
    m/sP*P*(n+2)*s;...
    (ddP - (n+1)*P)*c;...
    -(1/sP)*(m*dP + m*ctP*P)*s;...
    (dP*cP/sP - P*(m^2)/(sP^2) -P*(n+1))*c].*Bn;
    
    % Cnm expansion
    n = 2;
    m = 0;
    for pos = 2:Nc
        c = cos(m * lambda);
        s  = sin(m * lambda); 
        sP = sin(phic);
        cP = cos(phic);
        tP = sP/cP;
        ctP = 1/tP;
        P = double(subs(assocLegendre(n, m), phic));
        dP = double(subs(Diff_assocLegendre(n, m), phic));
        ddP = double(subs(Diff2_assocLegendre(n, m), phic));

        if(normalized)
            P = NormFactor(n, m) * P;
            dP = NormFactor(n, m) * dP;
            ddP = NormFactor(n, m) * ddP;
        end

        Bn = (GM/r^3) * (R/r)^n;

        GEOM(:, pos) = ...
            [(n + 1)*(n + 2)*P*c;...
            -dP*(n+2)*c;...
            m/sP*P*(n+2)*s;...
            (ddP - (n+1)*P)*c;...
            -(m/sP)*(dP + ctP*P)*s;...
            -(dP*cP/sP + P*(m^2)/(sP^2) +P*(n+1))*c].*Bn;

        if(m < n)
            m = m + 1;
        else
            m = 0;
            n = n + 1;
        end
    end

    % Snm expansion
    n = 2;
    m = 1;
    for pos = Nc+1:Nx
        c = cos(m * lambda);
        s  = sin(m * lambda); 
        sP = sin(phic);
        cP = cos(phic);
        tP = sP/cP;
        ctP = 1/tP;
        P = double(subs(assocLegendre(n, m), phic));
        dP = double(subs(Diff_assocLegendre(n, m), phic));
        ddP = double(subs(Diff2_assocLegendre(n, m), phic));

        if(normalized)
            P = NormFactor(n, m) * P;
            dP = NormFactor(n, m) * dP;
            ddP = NormFactor(n, m) * ddP;
        end

        Bn = (GM/r^3) * (R/r)^n;

        GEOM(:, pos) = ...
            [(n + 1)*(n + 2)*P*s;...
            -dP*(n+2)*s;...
            -m/sP*P*(n+2)*c;...
            (ddP - (n+1)*P)*s;...
            m/sP*(dP + ctP*P)*c;...
            -(dP*cP/sP + P*(m^2)/(sP^2) +P*(n+1))*s].*Bn;
  
        if(m < n)
            m = m + 1;
        else
            m = 1;
            n = n + 1;
        end
    end
end

%% FUNCTIONS
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
