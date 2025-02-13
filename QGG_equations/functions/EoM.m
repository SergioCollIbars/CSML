function [dx] = EoM(t, x, GM, OrbitObj)
    %%                          EoM FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 31/10/2022
    %
    %   Description: This function defines the equation of motion (EoM) for
    %   the orbital problem
    %
    %   Input:
    %       t: time vector
    %       x: state vector [r1, r2, r3, v1, v2, v3]'
    %       mu: gravity parameter
    %       OrbitObj: Orbit object
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%
    
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    [U] = potentialGradient_GM(GM, x1, x2, x3);
    
    % differential equations
       dx = [x(4);
          x(5);
          x(6);
          U(1);
          U(2);
          U(3)];

end