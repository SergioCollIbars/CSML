function [x,y,z] = fibonacci_sphere(N, R)
    % Generate N quasi-uniform points on a sphere of radius R
    %
    % Inputs:
    %   N : number of points
    %   R : sphere radius
    %
    % Outputs:
    %   x,y,z : Cartesian coordinates (Nx1)
    
    if nargin < 2
        R = 1;
    end
    
    k = (0:N-1)';
    golden_angle = pi * (3 - sqrt(5));
    
    theta = k * golden_angle;
    z = 1 - 2*(k + 0.5)/N;
    r_xy = sqrt(1 - z.^2);
    
    x = R * r_xy .* cos(theta);
    y = R * r_xy .* sin(theta);
    z = R * z;
end
