function [RTN_2_ECI] = RTN2ECI(r, v)
    % ------------------------------------------------------------------- %
    %                     RTN 2 ECI FRAME FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 02/09/2023

    % Description: Compute the radial, along-track, cross-track (RTN) 
    %   to Earth-centered inertial rotation matrix. If applied to a
    %   position vector in the RTN frame, it will transform that vector to
    %   into the equivalent position vector in the ECI frame.

    % Input:
    %   r: current position vector. ECI frame
    %   v: current velocity vector. ECI frame

    % Output:
    %   RTN_2_ECI: rotation matrix from RTN to ECI frame
    % ------------------------------------------------------------------- %

    n = cross(r, v);

    R = r / vecnorm(r);
    N = n / vecnorm(n);
    T = cross(N, R);

    RTN_2_ECI = [R, T, N];

end

