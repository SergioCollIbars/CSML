function [Hacc, Hgg] = potentialGradient_Cnm(n_max, r, Re, GM, ...
    NB, normalized)
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
    %       normalized: normalize coefficents. 1 / 0 
    %       NB: ACAF to ACI rotation matrix @ current time
    %
    %   Output: 
    %       Hgg: visibility matrix. pratials with respect to coeff for
    %       gradiometer measurements.
    %
    %       Hacc: visibility matrix. pratials with respect to coeff for
    %       acceleration measurements.
    % --------------------------------------------------------------------%
    
    % position magnitude. ECEF
    r_n = vecnorm(r);

    % ECEF coordinates
    x = r(1);
    y = r(2);
    z = r(3);

    % define visibility matrix
    Nxc = 0;
    for j = 1:n_max+1
        Nxc  = Nxc + j;
    end
    Nxc = Nxc - 2;
    Nxs = Nxc - n_max;
    Hgg = zeros(9, Nxc + Nxs);
    Hacc = zeros(3, Nxc + Nxs);

    % define b coefficient matrix
    b_real = zeros(n_max + 3, n_max + 3);
    b_imag = zeros(n_max + 3, n_max + 3);
        
    % Gradiometer spherical partials matrix components
    dxdx_dc = zeros(n_max + 1, n_max + 1) * NaN;
    dxdy_dc = zeros(n_max + 1, n_max + 1) * NaN;
    dxdz_dc = zeros(n_max + 1, n_max + 1) * NaN;

    dydy_dc = zeros(n_max + 1, n_max + 1) * NaN;
    dydz_dc = zeros(n_max + 1, n_max + 1) * NaN;

    dzdz_dc = zeros(n_max + 1, n_max + 1) * NaN;

    dxdx_ds = zeros(n_max + 1, n_max + 1) * NaN;
    dxdy_ds = zeros(n_max + 1, n_max + 1) * NaN;
    dxdz_ds = zeros(n_max + 1, n_max + 1) * NaN;

    dydy_ds = zeros(n_max + 1, n_max + 1) * NaN;
    dydz_ds = zeros(n_max + 1, n_max + 1) * NaN;

    dzdz_ds = zeros(n_max + 1, n_max + 1) * NaN;
    
    % Acceleration spherical partials matrix components.
    dx_dc = zeros(n_max + 1, n_max + 1) * NaN;
    dy_dc = zeros(n_max + 1, n_max + 1) * NaN;
    dz_dc = zeros(n_max + 1, n_max + 1) * NaN;

    dx_ds = zeros(n_max + 1, n_max + 1) * NaN;
    dy_ds = zeros(n_max + 1, n_max + 1) * NaN;
    dz_ds = zeros(n_max + 1, n_max + 1) * NaN;
    
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

            % compute g and s functions
            if(normalized == 1) % Normalized coefficients
                [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
                    getGS_normalized(n, m);
            elseif(normalized == 0) % Unnormalized coefficients
                [g1, g2, g3, s1, s2, s3, s4, s5, s6] = ...
                    getGS_unnormalized(n, m);
            end

            % acceleration spherical harmonics derivative components
            if(m == 0)
                dx_dc(N, M) = GM / (Re^2) * g1 * ...
                    (- b_real(N+1, M+1));

                dx_ds(N, M) = GM / (Re^2) * g1 * ...
                    (- b_imag(N+1, M+1));

                dy_dc(N, M) = GM / (Re^2) * g1 * ...
                    (- b_imag(N+1, M+1));

                dy_ds(N, M) = GM / (Re^2) * g1 * ...
                    (b_real(N+1, M+1));
            else
                dx_dc(N, M) = GM / (Re^2) * ...
                    (g1 * (-b_real(N+1, M+1)) + ...
                    g2 * (b_real(N+1, M-1)));

                dx_ds(N, M) = GM / (Re^2) * ...
                    (g1 * (-b_imag(N+1, M+1)) + ...
                    g2 * (b_imag(N+1, M-1)));
                
                dy_dc(N, M) = GM / (Re^2) * ...
                    (g1 * (-b_imag(N+1, M+1) ) + ...
                    g2 * ( -b_imag(N+1, M-1) ));

                dy_ds(N, M) = GM / (Re^2) * ...
                    (g1 * (b_real(N+1, M+1)) + ...
                    g2 * (b_real(N+1, M-1)));
            end
                dz_dc(N, M) = GM / (Re^2) * g3 * ...
                    (b_real(N+1, M));
                
                dz_ds(N, M) = GM / (Re^2) * g3 * ...
                    (b_imag(N+1, M));
    
            % gradiometer spherical harmonics derivative components
            if(m == 0)
                 dxdx_dc(N, M) = GM / (Re^3) * 0.5 * ((s1 * b_real(N + 2, 3) - ...
                     s2 * b_real(N + 2, 1)));

                 dxdx_ds(N, M) = 0;

                 dxdy_dc(N, M) = GM /(Re^3) * 0.5 * ...
                     (s1 * b_imag(N + 2, 3));

                 dxdy_ds(N, M) = 0;

            elseif( m == 1)
                dxdx_dc(N, M) = GM / (Re^3) * 0.25 * ...
                    (s1 * b_real(N + 2, 4) - 3 * s2 * ...
                    b_real(N + 2, 2));

                 dxdx_ds(N, M) = GM / (Re^3) * 0.25 * ...
                     ((s1* b_imag(N + 2, 4) - ...
                    s2 * b_imag(N + 2, 2)));

                 dxdy_dc(N, M) = GM / (Re^3) * 0.25 * ((s1 * b_imag(N + 2, 4) ...
                     - s2 * b_imag(N + 2, 2))); 
                
                 dxdy_ds(N, M) = GM / (Re^3) * 0.25 * ...
                     (( -(s1 * b_real(N + 2, 4) + s2 * b_real(N + 2, 2))));
                
            elseif(m > 1)
                dxdx_dc(N, M) = GM / (Re^3) * 0.25 * ...
                    (((s1 * b_real(N + 2, M + 2 ) - 2 * s2 * b_real(N + 2, M) + ...
                    s3 * b_real(N + 2, M - 2))));

                dxdx_ds(N, M) = GM / (Re^3) * 0.25 * ...
                    (((s1 *  b_imag(N + 2, M + 2) ...
                    - 2 * s2 * b_imag(N + 2, M) ...
                    + s3 * b_imag(N + 2, M - 2))));
                
                dxdy_dc(N, M) = GM / (Re^3) * 0.25 * ...
                     ((( s1 * b_imag(N + 2, M + 2) - s3 * b_imag(N + 2, M - 2))));

                dxdy_ds(N, M) = GM / (Re^3) * 0.25 * ...
                     ((-(s1 * b_real(N + 2, M + 2) - s3 * ...
                     b_real(N + 2, M - 2) )));
            end

            if(m == 0)
                dxdz_dc(N, M) = GM / (Re^3) * s4 * b_real(N + 2, 2);
                dxdz_ds(N, M) = 0;
 
                dydz_dc(N, M) = GM / (Re^3) * s4 * b_imag(N + 2, 2);
                dydz_ds(N, M) = 0;

            else
                dxdz_dc(N, M) = GM / (Re^3) * 0.5 * ...
                    ((( s4 *  b_real(N + 2, M + 1) -...
                    s5 * b_real(N + 2, M -1) )));

                dxdz_ds(N, M) = GM / (Re^3) * 0.5 * ...
                    ((( s4 *  ...
                    b_imag(N + 2, M + 1) - s5 * b_imag(N + 2, M - 1))));

                dydz_dc(N, M) = GM / (Re^3) * 0.5 * ...
                    (((s4 *...
                    b_imag(N + 2, M + 1) + s5 * b_imag(N +2 , M -1))));

                dydz_ds(N, M) = GM / (Re^3) * 0.5 * ...
                    ((-( s4 *  b_real(N + 2, M + 1) + ...
                    s5 * ...
                    b_real(N + 2, M - 1))));
                    
            end

            dzdz_dc(N, M) = GM / (Re^3) * (s6 * (b_real(N + 2, M)));

            dzdz_ds(N, M) = GM / (Re^3) * (s6 * (b_imag(N + 2, M)));

            % Use laplace eq. to compute last term
            dydy_dc(N, M) = - dxdx_dc(N, M) - dzdz_dc(N, M);
            dydy_ds(N, M) = - dxdx_ds(N, M) - dzdz_ds(N, M);
        end
    end
    % fill symetric values
    dydx_dc = dxdy_dc;
    dzdx_dc = dxdz_dc;
    dzdy_dc = dydz_dc;

    dydx_ds = dxdy_ds;
    dzdx_ds = dxdz_ds;
    dzdy_ds = dydz_ds;

    % rotate from ECEF to ECI frame
    for nl = 1:n_max+1
        for i = 1:nl
            a = [dxdx_dc(nl, i),dxdy_dc(nl, i),dxdz_dc(nl, i);...
                 dydx_dc(nl, i),dydy_dc(nl, i),dydz_dc(nl, i);...
                 dzdx_dc(nl, i),dzdy_dc(nl, i),dzdz_dc(nl, i)];
            a =  NB * a * NB';

            b = [dxdx_ds(nl, i),dxdy_ds(nl, i),dxdz_ds(nl, i);...
                 dydx_ds(nl, i),dydy_ds(nl, i),dydz_ds(nl, i);...
                 dzdx_ds(nl, i),dzdy_ds(nl, i),dzdz_ds(nl, i)];
            b =  NB * b * NB';
            
            dxdx_dc(nl, i) = a(1, 1);
            dxdy_dc(nl, i) = a(1, 2);
            dxdz_dc(nl, i) = a(1, 3);
            dydx_dc(nl, i) = a(2, 1);
            dydy_dc(nl, i) = a(2, 2);
            dydz_dc(nl, i) = a(2, 3);
            dzdx_dc(nl, i) = a(3, 1);
            dzdy_dc(nl, i) = a(3, 2);
            dzdz_dc(nl, i) = a(3, 3);

            dxdx_ds(nl, i) = b(1, 1);
            dxdy_ds(nl, i) = b(1, 2);
            dxdz_ds(nl, i) = b(1, 3);
            dydx_ds(nl, i) = b(2, 1);
            dydy_ds(nl, i) = b(2, 2);
            dydz_ds(nl, i) = b(2, 3);
            dzdx_ds(nl, i) = b(3, 1);
            dzdy_ds(nl, i) = b(3, 2);
            dzdz_ds(nl, i) = b(3, 3);
        end
    end

    % ensamble visibility matrix. C harmonics
    Hgg(:, 1) = [dxdx_dc(1, 1);dydx_dc(1, 1);dzdx_dc(1, 1);...
               dxdy_dc(1, 1);dydy_dc(1, 1);dzdy_dc(1, 1);...
               dxdz_dc(1, 1);dydz_dc(1, 1);dzdz_dc(1, 1)];

    Hacc(:, 1) = [dx_dc(1,1);dy_dc(1,1);dz_dc(1,1)];
    
    j = 2;
    nl = 3;
    while j < Nxc
        for i = 1:nl
            Hgg(1, j) = dxdx_dc(nl, i);
            Hgg(2, j) = dydx_dc(nl, i);
            Hgg(3, j) = dzdx_dc(nl, i);

            Hgg(4, j) = dxdy_dc(nl, i);
            Hgg(5, j) = dydy_dc(nl, i);
            Hgg(6, j) = dzdy_dc(nl, i);

            Hgg(7, j) = dxdz_dc(nl, i);
            Hgg(8, j) = dydz_dc(nl, i);
            Hgg(9, j) = dzdz_dc(nl, i);

            Hacc(1, j) = dx_dc(nl, i);
            Hacc(2, j) = dy_dc(nl, i);
            Hacc(3, j) = dz_dc(nl, i);
            Hacc(:, j) = NB * Hacc(:, j);

            j = j + 1;
        end
        nl = nl + 1;
    end

    j = Nxc + 1;
    nl = 3;
    while j < Nxc + Nxs
        for i = 1:nl
            if(i ~= 1)
                Hgg(1, j) = dxdx_ds(nl, i);
                Hgg(2, j) = dydx_ds(nl, i);
                Hgg(3, j) = dzdx_ds(nl, i);
    
                Hgg(4, j) = dxdy_ds(nl, i);
                Hgg(5, j) = dydy_ds(nl, i);
                Hgg(6, j) = dzdy_ds(nl, i);
    
                Hgg(7, j) = dxdz_ds(nl, i);
                Hgg(8, j) = dydz_ds(nl, i);
                Hgg(9, j) = dzdz_ds(nl, i);

                Hacc(1, j) = dx_ds(nl, i);
                Hacc(2, j) = dy_ds(nl, i);
                Hacc(3, j) = dz_ds(nl, i);
                Hacc(:, j) = NB * Hacc(:, j);

                j = j + 1;
            end
        end
        nl = nl + 1;
    end
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

