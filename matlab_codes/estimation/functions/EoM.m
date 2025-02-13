function [dx] = EoM(t, x, Cnm_t, Snm_t, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC)
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
    %       Cnm / Snm: SH coefficients for the asteroid (matrix form)
    %       n_max: maximum SH degree
    %       GM: gravitational parameter asteroid
    %       Re: asteroid reference radius
    %       normalized: normalized SH coefficientes: 1 yes / 0 no
    %       W0, W, RA, DEC: asteroid pole parameters
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%
    
     % position vector ACI
    r_ACI = [x(1); x(2); x(3)];

    % ACI to ACAF rotation matrix
    Wt = W0 + W * t;
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    r_ACAF = ACAF_ACI * r_ACI;

    [~, dU_ACAF, ddU_ACAF] = potentialGradient_nm(Cnm_t, Snm_t, n_max, ...
                                         r_ACAF, Re, GM, normalized);
    dU_ACI = ACAF_ACI' * dU_ACAF;

    ddU_ACI = ACAF_ACI' * ddU_ACAF * ACAF_ACI;
        
    J = [zeros(3, 3), eye(3, 3);ddU_ACI, zeros(3, 3)];

    PHI = reshape(x(7:end), [6, 6]);
    PHI_dot = J  * PHI;

    % differential equations
       dx = [x(4);
          x(5);
          x(6);
          dU_ACI(1);
          dU_ACI(2);
          dU_ACI(3);
          reshape(PHI_dot, [36, 1])];
end