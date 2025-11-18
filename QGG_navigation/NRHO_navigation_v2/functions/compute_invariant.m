function [I1,I2, I3] = compute_invariant(T)
    %   COMPUTE_INVARIANT 
    % Desfription: compute symmetric tensor invariant
    % Input: T: tensor value at time t
    % Output: Invariant: 1, 2, 3
    
    % get tensor components
    A11 = T(1, :);
    A12 = T(2, :);
    A13 = T(3, :);

    A22 = T(4, :);
    A23 = T(5, :);
    A33 = T(6, :);
    
    I1 = A11 + A22 + A33;
    I2 = A11.*A22 + A11.*A33 + A22.*A33 - A12.*A12 - A13.*A13 -A23.*A23;
    I3 = A11.*A22.*A33 + 2*A12.*A13.*A23 - A11.*A23.*A23 - A22.*A13.*A13 +...
        -A33.*A12.*A12;
end

