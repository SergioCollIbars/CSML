function [Hrot_ang] = compute_angularDyadPartials(omega, H_omega, H_omegaDot)
    % Description: Using the partials of the angular velocity and
    % acceleration w.r.t the Euler angles, compute the Dyad partials w.r.t
    % the Euler angles.
    
    % partials w.r.t angular acceleration
    H_omegaDot = [zeros(1, 3); - H_omegaDot(3, :); H_omegaDot(2, :); ...
        H_omegaDot(3, :); zeros(1, 3); - H_omegaDot(1, :); ...
        -H_omegaDot(2, :); H_omegaDot(1, :); zeros(1, 3)];

    % partials w.r.t angular velocity. For a (3-2-1) rotation
    dw11_dt = -2*omega(2)*H_omega(2, :) - 2*omega(3)*H_omega(3, :);
    dw12_dt = omega(2)*H_omega(1, :) + omega(1)*H_omega(2, :);
    dw13_dt = omega(3)*H_omega(1, :) + omega(1)*H_omega(3, :);
    dw22_dt = -2*omega(1)*H_omega(1, :) - 2*omega(3)*H_omega(3, :);
    dw23_dt = omega(3)*H_omega(2, :) + omega(2)*H_omega(3, :);
    dw33_dt = -2*omega(1)*H_omega(1, :) - 2*omega(2)*H_omega(2, :);

    H_omega = [dw11_dt; dw12_dt; dw13_dt; dw12_dt; dw22_dt; dw23_dt;...
        dw13_dt; dw23_dt; dw33_dt];
    Hrot_ang = H_omega + H_omegaDot;
end

