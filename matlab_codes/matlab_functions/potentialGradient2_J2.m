function [ddU] = potentialGradient2_J2(GM, J2, Re, x, y, z)
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
    r9 = r^9;

    dUxx = (-4*x^4 - 3*x^2 * (y^2 - 9*z^2) + y^4 - 3*y^2*z^2 -4*z^4) /r9;
    dUxy = - (5*x*y * (x^2 + y^2 - 6*z^2)) / r9;
    dUxz = - (5*x*z * (3*x^2 + 3*y^2 - 4*z^2)) / r9;

    dUyx = - (5*x*y * (x^2 + y^2 - 6*z^2)) / r9;
    dUyy = (x^4 - 3*x^2 * (y^2 + z^2) - 4*y^4 + 27*y^2*z^2 -4*z^4) /r9;
    dUyz = - (5*y*z * (3*x^2 + 3*y^2 - 4*z^2)) / r9;

    dUzx = - (5*x*z * (3*x^2 + 3*y^2 - 4*z^2)) / r9;
    dUzy = - (5*y*z * (3*x^2 + 3*y^2 - 4*z^2)) / r9;
    dUzz = (3*x^4 + 6*x^2 * (y^2 - 4*z^2) + 3*y^4 - 24*y^2*z^2 +8*z^4) /r9;

    ddU = -3/2 * GM * J2 * Re * Re * [dUxx, dUxy, dUxz;...
                                      dUyx, dUyy, dUyz;...
                                      dUzx, dUzy, dUzz];
end