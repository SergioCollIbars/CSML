function [J] = compute_jacobian(T, daSRP_dr, daSRP_deta)
    % Compute the Jacobian of the state. 
    % state: r, v, eta, bias
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % partials acceleration w.r.t bias
    da_db = zeros(3, 9);

    % partials bias rate w.r.t bias
    dbDot_db = zeros(9,9);

    J = [zeros(3, 3), eye(3,3), zeros(3, 1), zeros(3, 9); ...
        T+daSRP_dr, zeros(3, 3), daSRP_deta, da_db; ...
        zeros(1,16);...
        zeros(9, 7), dbDot_db];
end

