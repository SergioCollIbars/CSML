function [omega2M, OmegaM] = angularMatrix(omega, Omega)
    %%                      ANGULAR MATRIX
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 20/03/2023
    %
    %   Description: function to compute angular velocity and angular
    %   acceleration cross products operators given instant values.
    %
    %   Inputs: 
    %       omega: angular velocities @ current time
    %       Omega: angular acceleration @ current time
    %       
    %   Outputs:
    %       omega2M: dot product between cross product operator omega
    %       matrices.
    %       OmegaM: cross product operator for Omega values
    % --------------------------------------------------------------------%

    % cross product omega & Omega matrices
    OmegaM = [0, -Omega(3), Omega(2);...
              Omega(3), 0, -Omega(1);...
              -Omega(2), Omega(1), 0];

    omegaM = [0, -omega(3), omega(2);...
              omega(3), 0, -omega(1);...
              -omega(2), omega(1), 0];

    % dot product omegaM
    omega2M = omegaM * omegaM;
end