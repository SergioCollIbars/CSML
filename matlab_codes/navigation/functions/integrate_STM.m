function [dx] = integrate_STM(t, X, mu, x, y, z, r1, r2)
    PHI = reshape(X, [6, 6]);

    % compute SOGT
    T = gradmeas_rotFrame(mu,x, y, z, r1, r2);

    % compute Jacobian
    J = [zeros(3, 3), eye(3,3); T, [2,0,0;0,-2,0;0,0,0]];

    % time derivative of the STM
    PHI_dot = J * PHI;

    dx = reshape(PHI_dot, [36, 1]);
end

