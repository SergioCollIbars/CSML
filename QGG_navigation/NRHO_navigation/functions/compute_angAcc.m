function [angAcc_vec] = compute_angAcc(angVel_vec, dt)
    % Compute the angular accelration using finite differences.
    % dt : time step in seconds

    Nt         = length(angVel_vec(1, :));
    angAcc_vec = zeros(3, Nt);
    for k = 2:Nt-1
        angVel_next = angVel_vec(:, k+1);
        angVel_prev = angVel_vec(:, k-1);

        angAcc_vec(:, k) = (angVel_next - angVel_prev)./(2*dt); % [rad/s^2]
    end
end

