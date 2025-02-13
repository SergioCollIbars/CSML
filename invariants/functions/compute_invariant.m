function [I1,I2, I3] = compute_invariant(Gamma)
    %   COMPUTE_INVARIANT 
    % Desfription: compute symmetric tensor invariant
    % Input: Gamma: tensor value
    % Output: Invariant: 1, 2, 3
    
    % get tensor components
    A11 = Gamma(1, 1);
    A22 = Gamma(2, 2);
    A33 = Gamma(3, 3);

    A12 = Gamma(1, 2);
    A13 = Gamma(1, 3);
    A23 = Gamma(2, 3);
    
    I1 = A11 + A22 + A33;
    I2 = A11*A22 + A11*A33 + A22*A33 - A12*A12 - A13*A13 -A23*A23;
    I3 = A11*A22*A33 + 2*A12*A13*A23 - A11*A23*A23 - A22*A13*A13 +...
        -A33*A12*A12;
end

