function [v] = computeVel(r, t)
    % ------------------------------------------------------------------- %
    %                     COMPUTE VELOCITY FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 02/09/2023

    % Description: Compute the velocity in the inertial frame using finite
    % deferences

    % Input:
    %   r: position vector. ACI frame
    %   t: time vector

    % Output:
    %   v: velocity vector. ACI frame
    % ------------------------------------------------------------------- %

    % time steps
    Nt = length(t);
    
    % velocity vector
    v = zeros(3, Nt);
    
    v(:, 1) = (r(:, 2) - r(:, 1))./(t(2)-t(1));
    v(:, Nt) = (r(:, Nt) - r(:, Nt-1))./(t(Nt)-t(Nt-1)); 
    for k = 2:Nt - 1
        v(:, k) = (r(:, k+1) - r(:, k - 1))./(t(k+1) - t(k-1));
    end
end