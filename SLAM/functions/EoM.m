function [dx] = EoM(t, x, planetParams, poleParams, Cnm, Snm)
    %%                          EoM FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 18/11/2025
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
    %       augmented: augmented to the C and S coefficients.
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%
    
    % extract values
    GM          = planetParams(1); 
    Re          = planetParams(2);
    n_max       = planetParams(3);
    normalized  = planetParams(4);

    W           = poleParams(1);
    W0          = poleParams(2);
    RA          = poleParams(3);
    DEC         = poleParams(4);

    % position vector ACI
    r_ACI = [x(1); x(2); x(3)];

    % ACI to ACAF rotation matrix
    Wt = W0 + W * t;
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    r_ACAF = ACAF_ACI * r_ACI;

    [~, dU_ACAF, ddU_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                         r_ACAF, Re, GM, normalized);
    dU_ACI = ACAF_ACI' * dU_ACAF;

    ddU_ACI = ACAF_ACI' * ddU_ACAF * ACAF_ACI;

    [Nc, Ns, Ncs] = count_num_coeff(n_max);
    [Hacc, ~] = potentialGradient_Cnm(n_max, r_ACAF, Re, GM, ...
                        ACAF_ACI', normalized);
        
    Jxx = [zeros(3, 3), eye(3, 3), zeros(3,6);...
        ddU_ACI, zeros(3, 3),    zeros(3,6);...
        zeros(6, 12)];

    Jxc = [zeros(3, Nc), zeros(3, Ns);...
        Hacc;...
        zeros(6, Nc+Ns)];

    A = [Jxx, Jxc;...
        zeros(Ncs, 12), zeros(Ncs, Ncs)];

    Nx = 12 + Ncs;

    PHI = reshape(x(13:end), [Nx, Nx]);
    PHI_dot = A  * PHI;

    % differential equations
   dx = [x(4);
      x(5);
      x(6);
      dU_ACI(1);
      dU_ACI(2);
      dU_ACI(3);
      zeros(6, 1);
      reshape(PHI_dot, [Nx*Nx, 1])];
end