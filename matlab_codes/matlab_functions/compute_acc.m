function [ad] = compute_acc(G, M_E, omega, Omega, n, r, A1, A2, t, R)
    % Compute gravity gradient tensor

    % WARNING: This matrix is not constant. depends on the r position
    U_b = [2* G * M_E/ (r^3), 0, 0;
            0, -G * M_E / (r^3), 0;
            0, 0, -G * M_E / (r^3)];
   
    R3 = [cos(n) sin(n) 0;-sin(n) cos(n) 0;0 0 1];                 % B -> I
    U_i = R3 * U_b * R3';
    U = R' * U_i * R;

    % Construct angular acceleration matrixs
    Omega2 = [-omega(t, 3)^2-omega(t, 2)^2, omega(t, 1)*omega(t, 2), omega(t, 1)*omega(t, 3);
              omega(t, 1)*omega(t, 2), -omega(t, 3)^2-omega(t, 1)^2, omega(t, 2)*omega(t, 3);
              omega(t, 1)*omega(t, 3), omega(t, 2)*omega(t, 3), -omega(t, 1)^2-omega(t, 2)^2];
    
    OmegaDot = [0, -Omega(t, 3), Omega(t, 2);
                Omega(t, 3), 0, -Omega(t, 1);
                -Omega(t, 2), Omega(t, 1), 0];

    
    Delta_r1 =  A1';                  % Sensor 2 radius, Axis: {xs, ys, zs}
    Delta_r2 =  A2';                  % Sensor 2 radius, Axis: {xs, ys, zs}

    % Compute acceleration
    a1 = -(U - Omega2 - OmegaDot) * Delta_r1;
    a2 = -(U - Omega2 - OmegaDot) * Delta_r2;

    ad = (a1 - a2);

end