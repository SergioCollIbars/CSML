function [M] = compute_ExternalTorque(omega, I, dt)
    Nt = length(omega(1, :));
    M  = omega.*0;
    for k = 1:Nt-1
        angAcc_ap = (omega(:, k+1) - omega(:, k))./dt;
        A         = I * omega(:, k);
        M(:, k)   = I * angAcc_ap + cross(omega(:, k), A);
    end
end

