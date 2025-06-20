function [angAcc] = compute_angularVals_v2(angVel, I)
    % Description: given the Euler angles in a (3-2-1) rotation and the
    % Euler angles rate (velocity and acceleration) compute the angular
    % velocity and acceleration in the body frame. 
    % Compute the partials of the angular values with respect to the Euler
    % angles.
    Nt = length(angVel(1, :));
    angAcc = zeros(3, Nt);
    Iner = I;   % Inertia matrix
    for j = 1:Nt
        % angular acceleration
        t2 = Iner * angVel(:, j);
        angAcc(:, j) = -inv(Iner) * (cross(angVel(:, j), t2));
    end
end
