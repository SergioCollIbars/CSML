function [ad] = compute_acc_vec(omega, Omega, r, A1, A2, t, G, M_E)
    
    % create variables
    mu = G * M_E;
    omega = omega(:, t);
    Omega = Omega(:, t);

    % Acc sensor 1
    a1 = cross(omega, cross(omega, r)) + cross(omega, cross(omega, A1)) + ...
        cross(Omega, r) + cross(Omega, A1) + mu / vecnorm(r + A1)^3 * (r + A1);
    
    % Acc sensor 2
    a2 = cross(omega, cross(omega, r)) + cross(omega, cross(omega, A2)) + ...
        cross(Omega, r) + cross(Omega, A2) + mu / vecnorm(r + A2)^3 * (r + A2);
    
    % Relative acc
    ad = a1 - a2;

end