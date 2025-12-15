function [J] = compute_jacobian(T)
    % Compute the Jacobian of the state. 
    % state: r, v, bias
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J = [zeros(3, 3), eye(3,3); ...
        T, zeros(3, 3)];
end

