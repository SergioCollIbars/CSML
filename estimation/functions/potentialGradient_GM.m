function [dU] = potentialGradient_GM(GM, x, y, z)
    %%                    POTENTIAL GRADIENT GM
    % ------------------------------------------------------------------- %
    %
    % Author: Sergio Coll Ibars
    %
    % Date: 02/22/2023
    %
    % Description: Computes the potential gradient for GM.
    %
    % Inputs: 
    %   GM. Gravity parameter
    %   x coordinate
    %   y coordinate
    %   z coordinate
    % ------------------------------------------------------------------- %
    
    r = sqrt(x^2 + y^2 + z^2);
    r3 = r^3;

    dUx = x/r3;
    dUy = y/r3;
    dUz = z/r3;

    dU = - GM * [dUx; dUy; dUz];
end