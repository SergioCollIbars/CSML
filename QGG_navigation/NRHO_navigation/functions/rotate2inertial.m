function [r_I, v_I] = rotate2inertial(r_B, v_B, M, n)
    % M = wrapTo2Pi(M);
    N = length(M);
    r_I = zeros(3, N);
    v_I = zeros(3, N);
    W_tilde = [0,-n,0;n,0,0;0,0,0];
    W = [0;0;n];
    for j = 1:N
        % rotation matrix
        NB = [cos(M(j)), -sin(M(j)), 0;...
            sin(M(j)), cos(M(j)), 0;...
            0, 0, 1];
        NB_dot = - W_tilde * NB;

        % rotate to inertial
        r_I(:, j) = NB * r_B(:, j);
        v_I(:, j) = NB * v_B(:, j) + cross(W, r_I(:, j));
    end
end

