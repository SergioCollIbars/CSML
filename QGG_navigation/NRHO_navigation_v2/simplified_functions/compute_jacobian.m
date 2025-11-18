function [J] = compute_jacobian(T, daSRP_dr)
    % Compute the Jacobian of the state. 
    % state: r, v, eta, bias
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % partials acceleration w.r.t bias
    da_db = zeros(3, 6);

    % partials bias rate w.r.t bias
    dbDot_db = zeros(6,6);

    J = [zeros(3, 3), eye(3,3), zeros(3, 6); ...
        T+daSRP_dr, zeros(3, 3), da_db; ...
        zeros(6, 6), dbDot_db];
end

