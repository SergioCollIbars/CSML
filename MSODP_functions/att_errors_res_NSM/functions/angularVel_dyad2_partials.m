function [H_omega_stack] = angularVel_dyad2_partials(omega)
    H_omega = eye(3,3);

    dw11_dt = -2*omega(2)*H_omega(2, :) - 2*omega(3)*H_omega(3, :);
    dw12_dt = omega(2)*H_omega(1, :) + omega(1)*H_omega(2, :);
    dw13_dt = omega(3)*H_omega(1, :) + omega(1)*H_omega(3, :);
    dw22_dt = -2*omega(1)*H_omega(1, :) - 2*omega(3)*H_omega(3, :);
    dw23_dt = omega(3)*H_omega(2, :) + omega(2)*H_omega(3, :);
    dw33_dt = -2*omega(1)*H_omega(1, :) - 2*omega(2)*H_omega(2, :);

    H_omega_stack = [dw11_dt; dw12_dt; dw13_dt; dw12_dt; dw22_dt; dw23_dt;...
        dw13_dt; dw23_dt; dw33_dt];
end

