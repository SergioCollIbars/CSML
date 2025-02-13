function [dU] = potentialGradient_J2(GM, J2, Re, x, y, z)
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
    %   J2. 2nd harmonic coefficient
    %   Re. Planet radius
    %   x coordinate
    %   y coordinate
    %   z coordinate
    % ------------------------------------------------------------------- %
    
    r = sqrt(x^2 + y^2 + z^2);
    r7 = r^7;

    dUx = (x^2 + y^2 -4*z^2) * x / r7;
    dUy = (x^2 + y^2 -4*z^2) * y / r7;
    dUz = (3*x^2 +3*y^2 -2*z^2) * z / r7;

    dU = - 3/2 * GM * J2 * Re * Re * [dUx; dUy; dUz];
end