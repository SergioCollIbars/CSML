function [omegaDot] = compute_angAcc(omega, I, M)
    Nt = length(omega(1, :));
    omegaDot = omega.*0;
    for k = 1:Nt
        A = I * omega(:, k);
        omegaDot(:, k) = (I)\(M(:, k) - cross(omega(:, k), A));
    end
end

