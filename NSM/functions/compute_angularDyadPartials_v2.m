function [Hrot_omega_dyad, H_omegaDot_dyad, H_omega, H_omegaDot] = ...
    compute_angularDyadPartials_v2(omega, I)
    % Description: Using the partials of the angular velocity and
    % acceleration w.r.t the Euler angles, compute the Dyad partials w.r.t
    % the Euler angles. use the Euler eqns to compute angular acceleration
    %%%%%%%%%%%%%%%%%%%%%% PARTIALS EULER ANGLE %%%%%%%%%%%%%%%%%%%%%%%
    % partials w.r.t angular velocity. For a (3-2-1) rotation
    H_omega = eye(3,3);

    dw11_dt = -2*omega(2)*H_omega(2, :) - 2*omega(3)*H_omega(3, :);
    dw12_dt = omega(2)*H_omega(1, :) + omega(1)*H_omega(2, :);
    dw13_dt = omega(3)*H_omega(1, :) + omega(1)*H_omega(3, :);
    dw22_dt = -2*omega(1)*H_omega(1, :) - 2*omega(3)*H_omega(3, :);
    dw23_dt = omega(3)*H_omega(2, :) + omega(2)*H_omega(3, :);
    dw33_dt = -2*omega(1)*H_omega(1, :) - 2*omega(2)*H_omega(2, :);

    H_omega_stack = [dw11_dt; dw12_dt; dw13_dt; dw12_dt; dw22_dt; dw23_dt;...
        dw13_dt; dw23_dt; dw33_dt];

    Hrot_omega_dyad = H_omega_stack; % Dyad partial w.r.t angular velocity

    %%%%%%%%%%%%%%%%%%%%%% PARTIALS EULER RATE %%%%%%%%%%%%%%%%%%%%%%%
    % partials w.r.t angular velocity. For a (3-2-1) rotation

    Iner = I;    % inertia matrix
    
    [t1] = compute_skewMat(omega);
    [t2] = compute_skewMat(Iner * omega);
    H_omegaDot = -inv(Iner) * (t1 * Iner - t2);

    H_omegaDot_dyad = [zeros(1, 3); - H_omegaDot(3, :); H_omegaDot(2, :); ...
        H_omegaDot(3, :); zeros(1, 3); - H_omegaDot(1, :); ...
        -H_omegaDot(2, :); H_omegaDot(1, :); zeros(1, 3)]; % Dyad partials w.r.t angular velocity

end

function [skw] = compute_skewMat(vec)
    skw = [0, -vec(3), vec(2);...
           vec(3), 0, -vec(1);...
          -vec(2), vec(1), 0];
end