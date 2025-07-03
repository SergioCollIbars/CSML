function [U, dU, ddU] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                r, Re, GM, normalized)
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
    %       normalized: normalized option. 1 / 0
    %
    %   Output: 
    %       U: potential gradient in cartesian coordinates and ECEF frame
    %       dU: spacecraft acceleration in ECEF coordinates
    %       ddU: potential tensor
    %
    % --------------------------------------------------------------------%
    
    % position magnitude. ECEF
    r_n = vecnorm(r);

    % ECEF coordinates
    x = r(1);
    y = r(2);
    z = r(3);
    
    % define potential escalar
    U = 0;

    % define spacecraft (SC) acceleration. ECEF
    dU = zeros(3, 1);

    % A matrix components
    ddU_ddx = 0;
    ddU_ddz = 0;

    ddU_dxdy = 0;
    ddU_dxdz = 0;
    ddU_dydz = 0;
    
    % compute b values
    if(normalized == 1) % Normalized coefficients
        [b_real, b_imag] = getB_normalized(n_max, x, y, z, r_n, Re);
    elseif(normalized == 0) % Unnormalized coefficients
        [b_real, b_imag] = getB_unnormalized(n_max, x, y, z, r_n, Re);
    end

    % Compute acc, potential, A matrix
    for n = 0:n_max
        for m = 0:n
            % index delay
            N = n + 1;
            M = m + 1;

            % Current perturbed constants
            C = C_mat(N ,M);
            S = S_mat(N, M);

            % compute potential
            U = U + GM / Re * (b_real(N, M) * C + b_imag(N, M) * S);

            % compute g and s functions
            if(normalized == 1) % Normalized coefficients
                [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
                    getGS_normalized(n, m);
            elseif(normalized == 0) % Unnormalized coefficients
                [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
                    getGS_unnormalized(n, m);
            end
            
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
                 ddU_ddx = ddU_ddx +  GM / (Re^3) * 0.5 * (C_mat(N, 1) * ...
                     (s1 * b_real(N + 2, 3) - s2 * ...
                     b_real(N + 2, 1)));

                 ddU_dxdy = ddU_dxdy + GM /(Re^3) * 0.5 * ...
                     (C_mat(N, 1) * s1 * b_imag(N + 2, 3));

            elseif( m == 1)
                ddU_ddx = ddU_ddx +  GM / (Re^3) * 0.25 * ((C_mat(N, 2) *...
                    (s1 * b_real(N + 2, 4) - 3 * s2 * ...
                    b_real(N + 2, 2))) + (S_mat(N, 2) * (s1 * b_imag(N + 2, 4) - ...
                    s2 * b_imag(N + 2, 2)))); % WARNING. Change S 2nd index
                
                 ddU_dxdy = ddU_dxdy +  GM / (Re^3) * 0.25 * ...
                     (( - S_mat(N, 2) * (s1 * b_real(N + 2, 4) + s2 * b_real(N + 2, 2))) + (C_mat(N, 2) * ...
                     (s1 * b_imag(N + 2, 4) - s2 * b_imag(N + 2, 2)))); % WARNING. Change b_real 2nd index
                    
            elseif(m > 1)
                ddU_ddx = ddU_ddx +   GM / (Re^3) * 0.25 * ...
                    ((C * (s1 * b_real(N + 2, M + 2 ) - 2 * s2 * b_real(N + 2, M) + ...
                    s3 * ...
                    b_real(N + 2, M - 2))) + (S * (s1 * b_imag(N + 2, M + 2) ...
                    - 2 * s2 * b_imag(N + 2, M) ...
                    + s3 * b_imag(N + 2, M - 2))));
                
                 ddU_dxdy = ddU_dxdy + GM / (Re^3) * 0.25 * ...
                     (( -S * (s1 * b_real(N + 2, M + 2) - s3* ...
                     b_real(N + 2, M - 2) ) ) + (C * ...
                     (s1 * b_imag(N + 2, M + 2) - s3 * b_imag(N + 2, M - 2))));   
            end

            if(m == 0)
                ddU_dxdz = ddU_dxdz + GM / (Re^3) * s4 * C_mat(N, 1) * b_real(N + 2, 2);
                ddU_dydz = ddU_dydz + GM / (Re^3) * s4 * C_mat(N, 1) * b_imag(N + 2, 2);
            else
                ddU_dxdz = ddU_dxdz + GM / (Re^3) * 0.5 * ...
                    ((C *( s4 *  b_real(N + 2, M + 1) -...
                    s5 * ...
                    b_real(N + 2, M -1) ) ) + ( S * ( s4 *  ...
                    b_imag(N + 2, M + 1) - s5 * b_imag(N + 2, M - 1))));

                ddU_dydz = ddU_dydz + GM / (Re^3) * 0.5 * ...
                    ( ( - S * ( s4 *  b_real(N + 2, M + 1) + ...
                    s5 * ...
                    b_real(N + 2, M - 1))) + (C * (s4 *...
                    b_imag(N + 2, M + 1) + s5 * b_imag(N +2 , M -1))));
                    
            end

            ddU_ddz = ddU_ddz +  GM / (Re^3) * ( s6 * (C * b_real(N + 2, M) + S * b_imag(N + 2, M)));
        end
    end
    % Use laplace eq. to compute last term
    ddU_ddy = - ddU_ddx - ddU_ddz;

    % ensamble A matrix
    ddU = [ddU_ddx ddU_dxdy ddU_dxdz;...
         ddU_dxdy ddU_ddy ddU_dydz;...
         ddU_dxdz ddU_dydz ddU_ddz];

end

 %% FUNCTIONS
    function [b_real, b_imag] = getB_unnormalized(n_max, x, y, z, r_n, Re)
        % define b coefficient matrix
        b_real = zeros(n_max + 3, n_max + 3);
        b_imag = zeros(n_max + 3, n_max + 3);

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

    end
    % end function

    function [b_real, b_imag] = getB_normalized(n_max, x, y, z, r_n, Re)
        % define b coefficient matrix
        b_real = zeros(n_max + 3, n_max + 3);
        b_imag = zeros(n_max + 3, n_max + 3);
    
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
                        if n == 1
                            delta_1_n = 1.0;
                        else
                            delta_1_n = 0.0;
                        end
                        b_real(N, N) = sqrt((1.0 + delta_1_n) * ...
                            (2.0 * n + 1.0) / (2.0 * n)) * Re / (r_n) * ...
                            (x / r_n *b_real(N-1, N - 1) ...
                            - y / r_n * b_imag(N -1, N -1));
                        
                        b_imag(N, N) = sqrt((1.0 + delta_1_n) * ...
                            (2.0 * n + 1.0) / (2.0 * n)) * Re / (r_n) * ...
                            (y / r_n *b_real(N-1, N - 1) ...
                            + x / r_n * b_imag(N -1, N -1));
                    end
                else
                    
                    if (n >= 2)
                        b_real(N, M) =  sqrt((4.0 * n * n - 1.0) / (n * n - m * m)) * (Re * z) / (r_n^2) * ...
                            b_real(N-1, M) - sqrt((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m) / ((2.0 * n - 3.0) ...
                            * (n * n - m * m))) * (Re / r_n)^2 * ...
                            b_real(N-2, M);
        
                        b_imag(N, M) = sqrt((4.0 * n * n - 1.0) / (n * n - m * m)) * (Re * z) / (r_n^2) * ...
                            b_imag(N-1, M) - sqrt((2.0 * n + 1.0) * ((n - 1.0) * (n - 1.0) - m * m) / ((2.0 * n - 3.0) * ...
                            (n * n - m * m))) * (Re / r_n)^2 * ...
                            b_imag(N-2, M);
                    else
                        b_real(N, M) = sqrt((4.0 * n * n - 1.0) / (n * n - m * m)) * ...
                            (Re * z) / (r_n^2) * ...
                            b_real(N-1, M);
        
                        b_imag(N, M) = sqrt((4.0 * n * n - 1.0) / (n * n - m * m)) * ...
                            (Re * z) / (r_n^2) * ...
                            b_imag(N-1, M);
                    end
                end
            end
        end
    end
    % end function

    function [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
        getGS_unnormalized(n, m)
        % compute delta function
        if(m == 0) 
            delta = 1; 
        else 
            delta = 0; 
        end
            
        % compute g function
        g1 = 0.5 * (1 + delta);
        g2 = 0.5 * (n - m + 2) * (n - m + 1);
        g3 = - (n - m + 1);

        % compute s function
        s1 = 1;
        if(m == 0)
            s2 = (n + 2) * (n + 1);
        elseif(m == 1)
            s2 = (n + 1) * n;
        elseif(m > 1)
            s2 = (n - m + 2.0) * ...
                    (n - m + 1);
        end
        s3 = (n - m + 4) * (n - m + 3) * (n - m + 2) * (n - m + 1);
        if(m == 0)
            s4 = (n + 1); 
        else
            s4 = (n - m + 1);
        end
        s5 = (n - m + 3) * (n - m + 2) * (n - m + 1);
        s6 = (n - m + 2) * ...
                (n - m + 1);
    end
    % end function 

    function [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
        getGS_normalized(n, m)
        % compute delta function
        if m == 0
            delta_0_m = 1.0;
        else
            delta_0_m = 0.0;
        end % For if
        
        if m == 1
            delta_1_m = 1.0;
        else
            delta_1_m = 0.0;
        end % For if
        
        if m == 2
            delta_2_m = 1.0;
        else
            delta_2_m = 0.0;
        end
            
        % compute g function
        g1 = 0.5 * (1.0 + delta_0_m) * sqrt((2.0 - delta_0_m) * ...
            (2.0 * n + 1.0) * (n + m + 2.0) * ...
            (n + m + 1.0) / (2.0 * (2.0 * n + 3.0)));
        g2 = 0.5 * sqrt((2.0 - delta_0_m) * (2.0 * n + 1.0) * ...
            (n - m + 2.0) * (n - m + 1.0) / ((2.0 - delta_1_m) *...
            (2.0 * n + 3.0)));
        g3 = -sqrt((2.0 * n + 1.0) * ...
            (n + m + 1.0) * (n - m + 1.0) / (2.0 * n + 3.0));

        % compute s function
        s1 = sqrt(0.5 * (2.0 - delta_0_m) * (2.0 * n + 1.0) * ...
            (n + m + 4.0) * (n + m + 3.0) * (n + m + 2.0) * ...
            (n + m + 1.0) / (2.0 * n + 5.0));
        s2 = sqrt((2.0 * n + 1.0) * (n + m + 2.0) * (n + m + 1.0) * ...
            (n - m + 2.0) * (n - m + 1.0) / (2.0 * n + 5.0));
        s3 = sqrt((2.0 - delta_0_m) * (2.0 * n + 1.0) * (n - m + 4.0) * ...
            (n - m + 3.0) * (n - m + 2.0) * ...
            (n - m + 1.0) / ((2.0 - delta_2_m) * (2.0 * n + 5.0)));
        s4 = sqrt(0.5 * (2.0 - delta_0_m) * (2.0 * n + 1.0) * ...
            (n + m + 3.0) * (n + m + 2.0) * (n + m + 1.0) * ...
            (n - m + 1.0) / (2.0 * n + 5.0));
        s5 = sqrt((2.0 - delta_0_m) * (2.0 * n + 1.0) * ...
            (n + m + 1.0) * (n - m + 3.0) * (n - m + 2.0) * ...
            (n - m + 1.0) / ((2.0 - delta_1_m) * (2.0 * n + 5.0)));
        s6 = sqrt((2.0 * n + 1.0) * (n + m + 2.0) * (n + m + 1.0) * ...
            (n - m + 2.0) * (n - m + 1.0) / (2.0 * n + 5.0));
    end
    % end function 