function [U, dU, ddU] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                r, Re, GM, Normalized)
    %%                   COMPUTE ACCELERATION FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 01/18/2023
    %
    %   Description: This function computes the potential, potential 
    %   first and second gradient in the inertial frame.
    %
    %   Input:
    %       GM: gravity parameter
    %       Re: planet radious
    %       r: SC position vector ECEF (rotate)
    %       n_max: max n perturbation order
    %       C_mat: C coefficient matrix (non normalized)
    %       S_mat: S coefficient matrix (non normalized)
    %       Normalized: nomalized values for coefficients selector [1 / 0]
    %
    %   Output: 
    %       U: potential gradient in cartesian coordinates and ECEF frame
    %       dU: spacecraft acceleration in ECEF coordinates
    %       ddU: potential tensor ECEF coordinates
    %
    % --------------------------------------------------------------------%
    
    % position magnitude. ECEF
    r_n = vecnorm(r);

    % ECEF coordinates
    x = r(1);
    y = r(2);
    z = r(3);
    
    % define b coefficient matrix
    b_real = zeros(n_max + 3, n_max + 3);
    b_imag = zeros(n_max + 3, n_max + 3);
    
    % define potential escalar
    U = 0;

    % define spacecraft (SC) acceleration. ECEF
    dU = zeros(3, 1);

    % A matrix components
    ddU_ddx = 0;
    ddU_ddy = 0;
    ddU_ddz = 0;

    ddU_dxdy = 0;
    ddU_dxdz = 0;
    ddU_dydz = 0;
    
    % A_S matrix components
    ddx_dc = zeros(n_max + 1, n_max + 1) * NaN;
    ddy_dc = zeros(n_max + 1, n_max + 1) * NaN;
    ddz_dc = zeros(n_max + 1, n_max + 1) * NaN;

    ddx_ds = zeros(n_max + 1, n_max + 1) * NaN;
    ddy_ds = zeros(n_max + 1, n_max + 1) * NaN;
    ddz_ds = zeros(n_max + 1, n_max + 1) * NaN;
    
    % compute b values
    for m = 0:n_max + 2
        for n = m:n_max + 2

            % index delay
            N = n + 1;
            M = m + 1;

            % compute b coefficient
            if(m == n)
                if( m == 0)
                    b_real(N, M) = Re / r_n;
                    
                    b_imag(N, M) = 0;
                else
                    b_real(N, N) = (2*n -1) * Re / (r_n) * ...
                        (x / r_n *b_real(N-1, N - 1) ...
                        - y / r_n * b_imag(N -1, N -1));
                    
                    b_imag(N, N) = (2*n -1) * Re / (r_n) * (y / r_n *b_real(N-1, N - 1) ...
                        + x / r_n * b_imag(N -1, N -1));
                end
            else
                
                if (n >= 2)
                    b_real(N, M) = (2*n - 1) / (n - m) * (Re * z) / (r_n^2) * ...
                        b_real(N-1, M) - (n + m -1) / (n - m) * (Re / r_n)^2 * ...
                        b_real(N-2, M);
    
                    b_imag(N, M) = (2*n - 1) / (n - m) * (Re * z) / (r_n^2) * ...
                        b_imag(N-1, M) - (n + m -1) / (n - m) * (Re / r_n)^2 * ...
                        b_imag(N-2, M);
                else
                    b_real(N, M) = (2*n - 1) / (n - m) * (Re * z) / (r_n^2) * ...
                        b_real(N-1, M);
    
                    b_imag(N, M) = (2*n - 1) / (n - m) * (Re * z) / (r_n^2) * ...
                        b_imag(N-1, M);
                end
            end
        end
    end

    % Compute acc, potential, A matrix
    for n = 0:n_max
        if(Normalized)
            fac1 = factorial(n - 0);
            fac2 = factorial(n + 0);
            Nor1 = ((2 - 1)*(2*n + 1) * fac1 /fac2)^(0.5);
            if(n ~= 0)
                fac1 = factorial(n - 1);
                fac2 = factorial(n + 1);
                Nor2 = ((2 - 0)*(2*n + 1) * fac1 /fac2)^(0.5);
            else
                Nor2 = 1;
            end
        end
        for m = 0:n
            % index delay
            N = n + 1;
            M = m + 1;
            
            % compute delta function
            if(m == 0) 
                delta = 1; 
            else 
                delta = 0; 
            end

            % Current perturbed constants
            if(Normalized)
                fac1 = factorial(n - m);
                fac2 = factorial(n + m);
                Nor = ((2 - delta)*(2*n + 1) * fac1 /fac2)^(0.5);
            else
                 Nor = 1;
                 Nor1 = 1;
                 Nor2 = 1;
            end

            C = C_mat(N ,M) * Nor;
            S = S_mat(N, M) * Nor;

            C1 = C_mat(N, 1) * Nor1;
            C2 = C_mat(N, 2) * Nor2;

            S2 = S_mat(N, 2) * Nor2;

            % compute potential
            U = U + GM / Re * (b_real(N, M) * C + b_imag(N, M) * S);
            
            % compute g function
            g1 = 0.5 * (1 + delta);
            g2 = 0.5 * (n - m + 2) * (n - m + 1);
            g3 = - (n - m + 1);
            
            % compute SC acceleration ECEF
            if(m == 0)
                dU(1) = dU(1) +  GM / (Re^2) * g1 * ...
                    (- C * b_real(N+1, M+1) - ...
                    S * b_imag(N+1, M+1));

                dU(2) = dU(2) +  GM / (Re^2) * g1 * ...
                    (S * b_real(N+1, M+1) - ...
                    C * b_imag(N+1, M+1));
            else
                dU(1) = dU(1) + GM / (Re^2) * ...
                    (g1 * (-C * b_real(N+1, M+1) - S * b_imag(N+1, M+1) ) + ...
                    g2 * ( C * b_real(N+1, M-1) + S * b_imag(N+1, M-1) ));
                
                dU(2) = dU(2) + GM / (Re^2) * ...
                    (g1 * (S * b_real(N+1, M+1) - C * b_imag(N+1, M+1) ) + ...
                    g2 * ( S * b_real(N+1, M-1) - C * b_imag(N+1, M-1) ));
            end
                dU(3) = dU(3) + GM / (Re^2) * g3 * ...
                    ( C * b_real(N+1, M) + S * b_imag(N+1, M));
            
            

            % STM matrix components
            if(m == 0)
                 ddU_ddx = ddU_ddx +  GM / (Re^3) * 0.5 * (C1 * ...
                     (b_real(N + 2, 3) - (n + 2) * (n + 1) * ...
                     b_real(N + 2, 1)));

                 ddU_dxdy = ddU_dxdy + GM /(Re^3) * 0.5 * ...
                     (C1 * b_imag(N + 2, 3));

            elseif( m == 1)
                ddU_ddx = ddU_ddx +  GM / (Re^3) * 0.25 * ((C2 *...
                    (b_real(N + 2, 4) - 3 * (n + 1) * n * ...
                    b_real(N + 2, 2))) + (S2 * ( b_imag(N + 2, 4) - ...
                    (n + 1) * n * b_imag(N + 2, 2)))); % WARNING. Change S 2nd index
                
                 ddU_dxdy = ddU_dxdy +  GM / (Re^3) * 0.25 * ...
                     (( - S2 * ( b_real(N + 2, 4) + (n + 1) * ...
                     n * b_real(N + 2, 2))) + (C2 * ...
                     ( b_imag(N + 2, 4) - (n + 1) * n * b_imag(N + 2, 2)))); % WARNING. Change b_real 2nd index
                    
            elseif(m > 1)
                ddU_ddx = ddU_ddx +   GM / (Re^3) * 0.25 * ...
                    ((C * (b_real(N + 2, M + 2 ) - 2 * (n - m + 2.0) * ...
                    (n - m + 1) * b_real(N + 2, M) + (n - m + 4) * ...
                    (n - m + 3) * (n - m + 2) * (n - m + 1) * ...
                    b_real(N + 2, M - 2))) + (S * ( b_imag(N + 2, M + 2) ...
                    - 2 * (n - m + 2) * (n - m + 1) * b_imag(N + 2, M) ...
                    + (n - m + 4) * (n - m + 3) * (n - m + 2)...
                    * (n - m + 1) * b_imag(N + 2, M - 2))));
                
                 ddU_dxdy = ddU_dxdy + GM / (Re^3) * 0.25 * ...
                     (( -S * ( b_real(N + 2, M + 2) - (n - m + 4) * ...
                     (n - m + 3) * (n - m + 2) * (n - m + 1) * ...
                     b_real(N + 2, M - 2) ) ) + (C * ...
                     ( b_imag(N + 2, M + 2) - (n - m + 4) * (n - m + 3) *...
                     (n - m + 2) * (n - m + 1) * b_imag(N + 2, M - 2))));   
            end

            if(m == 0)
                ddU_dxdz = ddU_dxdz + GM / (Re^3) * (n + 1) * C1 * b_real(N + 2, 2);
                ddU_dydz = ddU_dydz + GM / (Re^3) * (n + 1) * C1 * b_imag(N + 2, 2);
            else
                ddU_dxdz = ddU_dxdz + GM / (Re^3) * 0.5 * ...
                    ((C *( (n - m + 1) *  b_real(N + 2, M + 1) -...
                    (n - m + 3) * (n - m + 2) * (n - m + 1) * ...
                    b_real(N + 2, M -1) ) ) + ( S * ( (n - m + 1) *  ...
                    b_imag(N + 2, M + 1) - (n - m + 3) * (n - m + 2) *...
                    (n - m + 1) * b_imag(N + 2, M - 1))));

                ddU_dydz = ddU_dydz + GM / (Re^3) * 0.5 * ...
                    ( ( - S * ( (n - m + 1) *  b_real(N + 2, M + 1) + ...
                    (n - m + 3) * (n - m + 2) * (n - m + 1) * ...
                    b_real(N + 2, M - 1))) + (C * ((n - m + 1) *...
                    b_imag(N + 2, M + 1) + (n - m + 3) * (n - m + 2) * ...
                    (n - m + 1) * b_imag(N +2 , M -1))));
                    
            end

            ddU_ddz = ddU_ddz +  GM / (Re^3) * ( (n - m + 2) * ...
                (n - m + 1) * (C * b_real(N + 2, M) + S * b_imag(N + 2, M)));

            % Spherical harmonics partials
            if( m ==0)
                ddx_dc(N, M) = GM / (Re^2) * g1 * (-b_real(N + 1, M + 1));
                ddy_dc(N, M) = GM / (Re^2) * g1 * (-b_imag(N + 1, M + 1));
            else
                ddx_dc(N, M) = GM / (Re^2) * (g1 * ...
                    ( - b_real(N +1, M + 1) ) + g2 * (b_real(N + 1, M -1)));
                
                ddy_dc(N, M) = GM / (Re^2) * (g1 * ...
                    ( - b_imag(N +1, M + 1) ) + g2 * (- b_imag(N + 1, M -1)));

                ddx_ds(N, M) = GM / (Re^2) * (g1 * ...
                    ( - b_imag(N +1, M + 1) ) + g2 * (b_imag(N + 1, M -1)));

                ddy_ds(N, M) = GM / (Re^2) * (g1 * ...
                    (b_real(N +1, M + 1) ) + g2 * (b_real(N + 1, M -1)));
            end
            ddz_dc(N, M) = GM / (Re^2) * g3 * b_real(N + 1, M);
            ddz_ds(N, M) = GM / (Re^2) * g3 * b_imag(N + 1, M);
        end
    end
    % Use laplace eq. to compute last term
    ddU_ddy = - ddU_ddx - ddU_ddz;

    % ensamble A matrix
    ddU = [ddU_ddx ddU_dxdy ddU_dxdz;...
         ddU_dxdy ddU_ddy ddU_dydz;...
         ddU_dxdz ddU_dydz ddU_ddz];

end