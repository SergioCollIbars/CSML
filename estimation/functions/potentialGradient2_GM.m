function [ddU] = potentialGradient2_GM(GM, x, y, z)
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
    r5 = r^5;

    dUxx = (1/r3) -3*(x^2)/r5;
    dUxy = -3 * (x*y/r5);
    dUxz = -3 * (x*z/r5);

    dUyx = -3*x*y/r5;
    dUyy = (1/r3) -3 * (y^2)/r5;
    dUyz = -3 * y*z/r5;

    dUzx = -3 * x*z/r5;
    dUzy = -3 * y*z/r5;
    dUzz = (1/r3) - 3 * (z^2)/r5;

    ddU = -GM *[dUxx, dUxy, dUxz;...
                dUyx, dUyy, dUyz;...
                dUzx, dUzy, dUzz];
end