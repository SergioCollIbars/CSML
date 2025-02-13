classdef AccClass
    %%                      ACCELERATION CLASS
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: Class to store all the information related with the
    %   acceleration module
    %
    % --------------------------------------------------------------------%
    properties
        % noise profile type
        noiseProfile
        
        % sensor gravity tensor
        accT
        
        % sensor STD. diagonal and out diagonal noise
        sigmaN_ii                       
        sigmaN_ij
        
        % Potential gradient. Body frame
        dUb
        
        % time vector
        t

        % saved data
        ri
        vi
        omega
        Omega
        Rtot
    end
end