function [Ax, Nx] = batchFilter(deltaY, Hi, Ri, Ax, Nx)
    %%                  BATCH FILTER FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 23/02/2023
    %
    %   Description: This function implements the batch filter and updates
    %   the normal equation at each time 
    %
    %   Input:
    %       deltaY:  O vector at current time.
    %       Hi: visibility matrix @ current time.
    %       Ri: measurement covariance at current time
    %       Ax: previous information matrix value
    %       Nx: previous N matrix value
    %
    %   Output:
    %       Ax:  new information matrix value
    %       Nx: new Nx matrix value
    % --------------------------------------------------------------------%
    
    Ri_inv = inv(Ri);

    Ax = Ax + (Hi' * Ri_inv * Hi);
    Nx = (Hi' * Ri_inv * deltaY) + Nx;

end

