function [Hk] = compute_posPartials(n_max, normalized, C_mat, S_mat, ...
                                                R, GM, r, ACAF_ACI, ACAF_B)
    %%                   COMPUTE Gradient Partials FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars and Evan Karen
    %
    %   Date: 02/25/2026
    %
    %   Description: This function computes the potential 
    %   third gradient in the inertial frame.
    %
    %   Input:
    %       GM: gravity parameter
    %       R: planet radious
    %       r: SC position vector ECEF (rotate)
    %       n_max: max n perturbation order
    %       C_mat: C coefficient matrix (non normalized)
    %       S_mat: S coefficient matrix (non normalized)
    %       normalized: normalized option. 1 / 0
    %       ACAF_ACI: Rotation matrix from ACI to ACAF
    %       ACAF_B: Rotation matrix from Body to ACAF
    %
    %   Output: 
    %       Hk: Sensitivity Matrix in ECEF Frame Where Hk (9x3) is:
    %        |G_xx|[(d/dx) (d/dy) (d/dz)]
    %        |G_xy|
    %        |G_xz|
    %        |G_yx|
    %        |G_yy|
    %        |G_yz|
    %        |G_zx|
    %        |G_zy|
    %        |G_zz|
    %
    % --------------------------------------------------------------------%
    
    ACI_ACAF = ACAF_ACI';
    B_ACAF = ACAF_B';

    r_pos = ACI_ACAF'*r;

    % position magnitude. ECEF
    r_n = vecnorm(r_pos);

    % ECEF coordinates
    x = r_pos(1);
    y = r_pos(2);
    z = r_pos(3);

    % Hk Matrix Components
    dddU_dddx = 0;
    dddU_ddxdy = 0;
    dddU_ddxdz = 0;
    dddU_dxddy = 0;
    dddU_dxdydz = 0;
    dddU_dddy = 0;
    dddU_ddydz = 0;
    dddU_dxddz = 0;
    dddU_dyddz = 0;
    dddU_dddz = 0;

    % compute b values
    if normalized
        [b_real, b_imag] = getB_normalized(n_max, x, y, z, r_n, R);
    else
        [b_real, b_imag] = getB_unnormalized(n_max, x, y, z, r_n, R);
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

            [T1,T2,T3,T4,T5,T6] = getT(n,m);

            if (normalized)
                [N1,N2,N3,N4,N5,N6,N7] = getN(n,m);
            else
                N1 = 1;
                N2 = 1;
                N3 = 1;
                N4 = 1;
                N5 = 1;
                N6 = 1;
                N7 = 1;
            end

            if m > 2
                dddU_dddx = dddU_dddx + (GM / (8 * R^4)) * (...
                    C*(-N1*b_real(N+3,M+3) + 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) + T6*N7*b_real(N+3,M-3))...
                    + S*(-N1*b_imag(N+3,M+3) + 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) + T6*N7*b_imag(N+3,M-3)));

                dddU_ddxdy = dddU_ddxdy + (GM / (8 * R^4)) * (...
                    S*(T6*N7*b_real(N+3,M-3) - T4*N5*b_real(N+3,M-1) - T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    - C*(T6*N7*b_imag(N+3,M-3) - T4*N5*b_imag(N+3,M-1) - T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_ddxdz = dddU_ddxdz + (GM / (4 * R^4)) * (...
                    C*(-T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) - T5*N6*b_real(N+3,M-2))...
                    + S*(-T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) - T5*N6*b_imag(N+3,M-2)));

                dddU_dxddy = dddU_dxddy + (GM / (8 * R^4)) * (...
                    C*(-T6*N7*b_real(N+3,M-3) - T4*N5*b_real(N+3,M-1) + T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    + S*(-T6*N7*b_imag(N+3,M-3) - T4*N5*b_imag(N+3,M-1) + T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_dxdydz = dddU_dxdydz + (GM / (4 * R^4)) * (...
                    S*(T1*N2*b_real(N+3,M+2) - T5*N6*b_real(N+3,M-2))...
                    - C*(T1*N2*b_imag(N+3,M+2) - T5*N6*b_imag(N+3,M-2)));

                dddU_dddy = dddU_dddy + (GM / (8 * R^4)) * (...
                    S*(-N1*b_real(N+3,M+3) - 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) - T6*N7*b_real(N+3,M-3))...
                    - C*(-N1*b_imag(N+3,M+3) - 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) - T6*N7*b_imag(N+3,M-3)));

                dddU_ddydz = dddU_ddydz + (GM / (4 * R^4)) * (...
                    C*(T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) + T5*N6*b_real(N+3,M-2))...
                    + S*(T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) + T5*N6*b_imag(N+3,M-2)));

                dddU_dxddz = dddU_dxddz + (GM / (2 * R^4)) * (...
                    C*(-T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    + S*(-T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));

                dddU_dyddz = dddU_dyddz + (GM / (2 * R^4)) * (...
                    S*(T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    - C*(T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));
            elseif m == 2
                dddU_dddx = dddU_dddx + (GM / (8 * R^4)) * (...
                    C*(-N1*b_real(N+3,M+3) + 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) - (n+2)*(n+1)*(n)*(n-1)*N5*b_real(N+3,M-1))...
                    + S*(-N1*b_imag(N+3,M+3) + 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) + (n+2)*(n+1)*(n)*(n-1)*N5*b_imag(N+3,M-1)));

                dddU_ddxdy = dddU_ddxdy + (GM / (8 * R^4)) * (...
                    S*(-(n+2)*(n+1)*(n)*(n-1)*N5*b_real(N+3,M-1) - T4*N5*b_real(N+3,M-1) - T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    - C*((n+2)*(n+1)*(n)*(n-1)*N5*b_imag(N+3,M-1) - T4*N5*b_imag(N+3,M-1) - T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_ddxdz = dddU_ddxdz + (GM / (4 * R^4)) * (...
                    C*(-T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) - T5*N6*b_real(N+3,M-2))...
                    + S*(-T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) - T5*N6*b_imag(N+3,M-2)));

                dddU_dxddy = dddU_dxddy + (GM / (8 * R^4)) * (...
                    C*((n+2)*(n+1)*(n)*(n-1)*N5*b_real(N+3,M-1) - T4*N5*b_real(N+3,M-1) + T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    + S*(-(n+2)*(n+1)*(n)*(n-1)*N5*b_imag(N+3,M-1) - T4*N5*b_imag(N+3,M-1) + T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_dxdydz = dddU_dxdydz + (GM / (4 * R^4)) * (...
                    S*(T1*N2*b_real(N+3,M+2) - T5*N6*b_real(N+3,M-2))...
                    - C*(T1*N2*b_imag(N+3,M+2) - T5*N6*b_imag(N+3,M-2)));

                dddU_dddy = dddU_dddy + (GM / (8 * R^4)) * (...
                    S*(-N1*b_real(N+3,M+3) - 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) + (n+2)*(n+1)*(n)*(n-1)*N5*b_real(N+3,M-1))...
                    - C*(-N1*b_imag(N+3,M+3) - 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) - (n+2)*(n+1)*(n)*(n-1)*N5*b_imag(N+3,M-1)));

                 dddU_ddydz = dddU_ddydz + (GM / (4 * R^4)) * (...
                    C*(T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) + T5*N6*b_real(N+3,M-2))...
                    + S*(T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) + T5*N6*b_imag(N+3,M-2)));

                 dddU_dxddz = dddU_dxddz + (GM / (2 * R^4)) * (...
                    C*(-T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    + S*(-T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));

                 dddU_dyddz = dddU_dyddz + (GM / (2 * R^4)) * (...
                    S*(T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    - C*(T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));
            elseif m == 1
                dddU_dddx = dddU_dddx + (GM / (8 * R^4)) * (...
                    C*(-N1*b_real(N+3,M+3) + 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) + (n+1)*(n)*N3*b_real(N+3,M+1))...
                    + S*(-N1*b_imag(N+3,M+3) + 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) - (n+1)*(n)*N3*b_imag(N+3,M+1)));

                dddU_ddxdy = dddU_ddxdy + (GM / (8 * R^4)) * (...
                    S*((n+1)*(n)*N3*b_real(N+3,M+1) - T4*N5*b_real(N+3,M-1) - T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    - C*(-(n+1)*(n)*N3*b_imag(N+3,M+1) - T4*N5*b_imag(N+3,M-1) - T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_ddxdz = dddU_ddxdz + (GM / (4 * R^4)) * (...
                    C*(-T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) + (n+2)*(n+1)*(n)*N4*b_real(N+3,M))...
                    + S*(-T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) - (n+2)*(n+1)*(n)*N4*b_imag(N+3,M)));

                dddU_dxddy = dddU_dxddy + (GM / (8 * R^4)) * (...
                    C*(-(n+1)*(n)*N3*b_real(N+3,M+1) - T4*N5*b_real(N+3,M-1) + T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    + S*((n+1)*(n)*N3*b_imag(N+3,M+1) - T4*N5*b_imag(N+3,M-1) + T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_dxdydz = dddU_dxdydz + (GM / (4 * R^4)) * (...
                    S*(T1*N2*b_real(N+3,M+2) + (n+2)*(n+1)*(n)*N4*b_real(N+3,M))...
                    - C*(T1*N2*b_imag(N+3,M+2) - (n+2)*(n+1)*(n)*N4*b_imag(N+3,M)));

                dddU_dddy = dddU_dddy + (GM / (8 * R^4)) * (...
                    S*(-N1*b_real(N+3,M+3) - 3*T2*N3*b_real(N+3,M+1) - 3*T4*N5*b_real(N+3,M-1) - (n+1)*(n)*N3*b_real(N+3,M+1))...
                    - C*(-N1*b_imag(N+3,M+3) - 3*T2*N3*b_imag(N+3,M+1) - 3*T4*N5*b_imag(N+3,M-1) + (n+1)*(n)*N3*b_imag(N+3,M+1)));

                 dddU_ddydz = dddU_ddydz + (GM / (4 * R^4)) * (...
                    C*(T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) - (n+2)*(n+1)*(n)*N4*b_real(N+3,M))...
                    + S*(T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) + (n+2)*(n+1)*(n)*N4*b_imag(N+3,M)));

                 dddU_dxddz = dddU_dxddz + (GM / (2 * R^4)) * (...
                    C*(-T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    + S*(-T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));

                 dddU_dyddz = dddU_dyddz + (GM / (2 * R^4)) * (...
                    S*(T2*N3*b_real(N+3,M+1) + T4*N5*b_real(N+3,M-1))...
                    - C*(T2*N3*b_imag(N+3,M+1) + T4*N5*b_imag(N+3,M-1)));
            elseif m == 0
                dddU_dddx = dddU_dddx + (GM / (8 * R^4)) * (...
                    C*(-N1*b_real(N+3,M+3) + 3*T2*N3*b_real(N+3,M+1) + 3*(n+2)*(n+1)*N3*b_real(N+3,M+1) - N1*b_real(N+3,M+3))...
                    + S*(-N1*b_imag(N+3,M+3) + 3*T2*N3*b_imag(N+3,M+1) - 3*(n+2)*(n+1)*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_ddxdy = dddU_ddxdy + (GM / (8 * R^4)) * (...
                    S*(-N1*b_real(N+3,M+3) + (n+2)*(n+1)*N3*b_real(N+3,M+1) - T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    - C*(N1*b_imag(N+3,M+3) - (n+2)*(n+1)*N3*b_imag(N+3,M+1) - T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_ddxdz = dddU_ddxdz + (GM / (4 * R^4)) * (...
                    C*(-T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) - (n+1)*N2*b_real(N+3,M+2))...
                    + S*(-T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) + (n+1)*N2*b_imag(N+3,M+2)));

                dddU_dxddy = dddU_dxddy + (GM / (8 * R^4)) * (...
                    C*(N1*b_real(N+3,M+3) + (n+2)*(n+1)*N3*b_real(N+3,M+1) + T2*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    + S*(-N1*b_imag(N+3,M+3) - (n+2)*(n+1)*N3*b_imag(N+3,M+1) + T2*N3*b_imag(N+3,M+1) + N1*b_imag(N+3,M+3)));

                dddU_dxdydz = dddU_dxdydz + (GM / (4 * R^4)) * (...
                    S*(T1*N2*b_real(N+3,M+2) - (n+1)*N2*b_real(N+3,M+2))...
                    - C*(T1*N2*b_imag(N+3,M+2) + (n+1)*N2*b_imag(N+3,M+2)));

                dddU_dddy = dddU_dddy + (GM / (8 * R^4)) * (...
                    S*(-N1*b_real(N+3,M+3) - 3*T2*N3*b_real(N+3,M+1) + 3*(n+2)*(n+1)*N3*b_real(N+3,M+1) + N1*b_real(N+3,M+3))...
                    - C*(-N1*b_imag(N+3,M+3) - 3*T2*N3*b_imag(N+3,M+1) - 3*(n+2)*(n+1)*N3*b_imag(N+3,M+1) - N1*b_imag(N+3,M+3)));

                dddU_ddydz = dddU_ddydz + (GM / (4 * R^4)) * (...
                    C*(T1*N2*b_real(N+3,M+2) + 2*T3*N4*b_real(N+3,M) + (n+1)*N2*b_real(N+3,M+2))...
                    + S*(T1*N2*b_imag(N+3,M+2) + 2*T3*N4*b_imag(N+3,M) - (n+1)*N2*b_imag(N+3,M+2)));

                dddU_dxddz = dddU_dxddz + (GM / (2 * R^4)) * (...
                    C*(-T2*N3*b_real(N+3,M+1) - (n+2)*(n+1)*N3*b_real(N+3,M+1))...
                    + S*(-T2*N3*b_imag(N+3,M+1) + (n+2)*(n+1)*N3*b_imag(N+3,M+1)));

                dddU_dyddz = dddU_dyddz + (GM / (2 * R^4)) * (...
                    S*(T2*N3*b_real(N+3,M+1) - (n+2)*(n+1)*N3*b_real(N+3,M+1))...
                    - C*(T2*N3*b_imag(N+3,M+1) + (n+2)*(n+1)*N3*b_imag(N+3,M+1)));
            end
            dddU_dddz = dddU_dddz - (GM / (R^4))*T3*(C*N4*b_real(N+3,M) + S*N4*b_imag(N+3,M));
        end
    end

    T_inertial(:,:,1) = [dddU_dddx,dddU_ddxdy,dddU_ddxdz;...
        dddU_ddxdy,dddU_dxddy,dddU_dxdydz;...
        dddU_ddxdz,dddU_dxdydz,dddU_dxddz];
    T_inertial(:,:,2) = [dddU_ddxdy,dddU_dxddy,dddU_dxdydz;...
        dddU_dxddy,dddU_dddy,dddU_ddydz;...
        dddU_dxdydz,dddU_ddydz,dddU_dyddz];
    T_inertial(:,:,3) = [dddU_ddxdz,dddU_dxdydz,dddU_dxddz;...
        dddU_dxdydz,dddU_ddydz,dddU_dyddz;...
        dddU_dxddz,dddU_dyddz,dddU_dddz];
    % T_body = T_inertial;

    % % Chat-GPT Rotation Method 1 (Slow) 
    T_body = zeros(3,3,3);
    for a = 1:3
        for b = 1:3
            for c = 1:3
                for i = 1:3
                    for j = 1:3
                        for k = 1:3
                            T_body(a,b,c) = T_body(a,b,c) + ...
                                B_ACAF(a,i)*B_ACAF(b,j)*B_ACAF(c,k)*T_inertial(i,j,k);
                        end
                    end
                end
            end
        end
    end

    % % Chat_GPT Method 3
    % T_body = tensorprod(B_ACI, tensorprod(B_ACI, tensorprod(B_ACI, T_inertial, 2, 1), 2, 2), 2, 3);

    Hk = [T_body(1,1,1),T_body(1,1,2),T_body(1,1,3);...
        T_body(1,2,1),T_body(1,2,2),T_body(1,2,3);...
        T_body(1,3,1),T_body(1,3,2),T_body(1,3,3);...
        T_body(2,1,1),T_body(2,1,2),T_body(2,1,3);...
        T_body(2,2,1),T_body(2,2,2),T_body(2,2,3);...
        T_body(2,3,1),T_body(2,3,2),T_body(2,3,3);...
        T_body(3,1,1),T_body(3,1,2),T_body(3,1,3);...
        T_body(3,2,1),T_body(3,2,2),T_body(3,2,3);...
        T_body(3,3,1),T_body(3,3,2),T_body(3,3,3)];

end

 %% FUNCTIONS

function [b_real, b_imag] = getB_unnormalized(n_max, x, y, z, r_n, Re)
    % define b coefficient matrix
    b_real = zeros(n_max + 4, n_max + 4);
    b_imag = zeros(n_max + 4, n_max + 4);

    % compute b values
    for m = 0:n_max + 3
        for n = m:n_max + 3

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
    b_real = zeros(n_max + 4, n_max + 4);
    b_imag = zeros(n_max + 4, n_max + 4);

    % compute b values
    for m = 0:n_max + 3
        for n = m:n_max + 3

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

function [T1,T2,T3,T4,T5,T6] = ...
        getT(n, m)
    T1 = (n-m+1);
    T2 = (n-m+2)*(n-m+1);
    T3 = (n-m+3)*(n-m+2)*(n-m+1);
    T4 = (n-m+4)*(n-m+3)*(n-m+2)*(n-m+1);
    T5 = (n-m+5)*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1);
    T6 = (n-m+6)*(n-m+5)*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1);
end

function [N1,N2,N3,N4,N5,N6,N7] = getN(n,m)
    delta0 = 0;
    delta1 = 0;
    delta2 = 0;
    delta3 = 0;
    if m==0
        delta0 = 1;
    elseif m==1
        delta1 = 1;
    elseif m==2
        delta2 = 1;
    elseif m==3
        delta3 = 1;
    end

    N1 = sqrt(((2-delta0)*(2*n + 1)*(n+m+6)*(n+m+5)*(n+m+4)*(n+m+3)*(n+m+2)*(n+m+1))/...
        (2*(2*n + 7)));

    N2 = sqrt(((2-delta0)*(2*n + 1)*(n+m+5)*(n+m+4)*(n+m+3)*(n+m+2)*(n+m+1))/...
        (2*(2*n + 7)*(n-m+1)));

    N3 = sqrt(((2-delta0)*(2*n + 1)*(n+m+4)*(n+m+3)*(n+m+2)*(n+m+1))/...
        (2*(2*n + 7)*(n-m+2)*(n-m+1)));

    N4 = sqrt(((2*n + 1)*(n+m+3)*(n+m+2)*(n+m+1))/...
        ((2*n + 7)*(n-m+3)*(n-m+2)*(n-m+1)));

    N5 = sqrt(((2-delta0)*(2*n + 1)*(n+m+2)*(n+m+1))/...
        ((2-delta1)*(2*n + 7)*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)));

    N6 = sqrt(((2-delta0)*(2*n + 1)*(n+m+1))/...
        ((2-delta2)*(2*n + 7)*(n-m+5)*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)));
    
    N7 = sqrt(((2-delta0)*(2*n + 1))/...
        ((2-delta3)*(2*n + 7)*(n-m+6)*(n-m+5)*(n-m+4)*(n-m+3)*(n-m+2)*(n-m+1)));
end
% end function
