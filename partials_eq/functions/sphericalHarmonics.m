function [U, dU, ddU] = sphericalHarmonics(GM, Re, n_max, C_mat, S_mat, rk, ...
                        option)
    %%                   SPHERICAL HARMONICS FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 03/15/2023
    %
    %   Description: This function computes the analytical expresion 
    %       for the spherical harmonics evaluated at the current point
    %
    %   Input:
    %       GM: gravity parameter
    %       Re: planet radious
    %       n_max: max n perturbation order
    %       C_mat: C coefficient matrix (non normalized)
    %       S_mat: S coefficient matrix (non normalized)
    %       rk: position norm to evaluate expression. ECEF coords
    %       option: option vector [output type, coordiantes, normalized]
    %
    %   Output: 
    %       U: potential gradient in cartesian coordinates and ECEF frame
    %       dU: spacecraft acceleration in cartesian coords, ECEF frame
    %       ddU: potential tensor. cartesian coords, ECEF frame
    %
    % --------------------------------------------------------------------%
    
    % options
    outType = option(1);
    coords = option(2);
    normalized = option(3);

    % define analytical variables. Spherical coordinates
    syms r phi lambda;
    syms x y z;

    % symbolic relation. Spherical 2 cartesian
    rc = sqrt(x^2 + y^2 + z^2);
    phic = atan2(z, sqrt(x^2 + y^2));
    lambdac = atan2(y, x);
    
    % compute b constant
    b = GM/r;

    % compute gravity potential
    U = 0;
    for n = 0:n_max
        for m = 0:n
            a = C_mat(n+1, m+1) * cos(m * lambda) + ...
                S_mat(n+1, m+1) * sin(m * lambda);

            % associated legendre polinomial
            Pnm = assocLegendre(n, m);
            if(normalized == "1")
                [N] = NormFactor(n, m);
                Pnm = N * Pnm;
            elseif(normalized == "0")
                Pnm = 1 * Pnm;
            end

            U = U + (Re^n)/(r^n) * Pnm * a; 
        end
    end
    U = b * U;
    
    if(coords == "cartesian")
        % convert from spherical 2 cartesian
        U = subs(U, {r, phi, lambda}, {rc, phic, lambdac});

        % compute potential gradient
        dU = [diff(U, x); diff(U, y); diff(U, z)];
    
        % compute potential tensor
        ddU = [diff(dU(1), x), diff(dU(1), y), diff(dU(1), z);...
               diff(dU(2), x), diff(dU(2), y), diff(dU(2), z);...
               diff(dU(3), x), diff(dU(3), y), diff(dU(3), z)];
    elseif(coords == "spherical")
        % compute potential gradient
        dU = [diff(U, r); diff(U, phi); diff(U, lambda)];
    
        % compute potential tensor
        ddU = [diff(dU(1), r), diff(dU(1), phi), diff(dU(1), lambda);...
               diff(dU(2), r), diff(dU(2), phi), diff(dU(2), lambda);...
               diff(dU(3), r), diff(dU(3), phi), diff(dU(3), lambda)];
    end

    if(outType == "numerical")
        if(coords == "cartesian")
            % evaluate at current point
            U = double(subs(U, {x, y, z}, {rk(1), rk(2), rk(3)}));
        
            dU = double(subs(dU, {x, y, z}, {rk(1), rk(2), rk(3)}));
        
            ddU = double(subs(ddU, {x, y, z}, {rk(1), rk(2), rk(3)}));
        elseif(coords == "spherical")
            rc = sqrt(rk(1)^2 + rk(2)^2 + rk(3)^2);
            phic = atan2(z, sqrt(rk(1)^2 + rk(2)^2));
            lambdac = atan2(rk(2), rk(1));

             % evaluate at current point
            U = double(subs(U, {r, phi, lambda}, {rc, phic, lambdac}));
        
            dU = double(subs(dU, {r, phi, lambda}, {rc, phic, lambdac}));
        
            ddU = double(subs(ddU, {r, phi, lambda}, {rc, phic, lambdac}));
        end
    end
end

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