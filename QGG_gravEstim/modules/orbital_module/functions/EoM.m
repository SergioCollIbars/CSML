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
    %       GM: gravity parameter
    %       OrbitObj: Orbit object
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%

    % Get orbit data
    Re = OrbitObj.OrbitData(10);
    C_mat = OrbitObj.C;
    S_mat = OrbitObj.S;
    n_max = OrbitObj.OrbitData(9);
    Normalized = OrbitObj.OrbitData(13);

    W = OrbitObj.W;
    W0 = OrbitObj.W0;
    RA = OrbitObj.RA;
    DEC = OrbitObj.DEC;

    % ACI position vector
    r_ACI = [x(1); x(2); x(3)];
    v_ACI = [x(4); x(5); x(6)];

    % Inertial to ACAF rotation 
    Wt = W * t + W0;
    ACAF_N = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    % ACAF position vector
    r_ACAF = ACAF_N * r_ACI;
    
    % compute S/C acc
    [~, dU, ~] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                r_ACAF, Re, GM, ...
                                                Normalized);
    % rotate from ECEF 2 ECI
    dU = ACAF_N' * dU;
    
    % compute planet keplerian motion around Sun
    [rs_ACI, vs_ACI] = Kepler_motion(OrbitObj.alphaK, OrbitObj.GMs, t);
    H = cross(rs_ACI, vs_ACI);
    H = H / vecnorm(H);
    
    D = rs_ACI./vecnorm(rs_ACI);
    F = cross(H, D);

    SAR_N = [D';F';H'];

    % SRP model. Inertial frame
    if(OrbitObj.OrbitData(11) == 1)
        rc = rs_ACI;
        rs = rc + r_ACI;
        acc_SRP = SRP(OrbitObj, rs);

        dU = dU + acc_SRP;
    end

    % 3BP model. Inertial frame
    if(OrbitObj.OrbitData(12) == 1)
        semiMajor = OrbitObj.alphaK(1);
        F = sqrt(OrbitObj.GMs / semiMajor^3);
        r_SAR = SAR_N * r_ACI;
        BP_x = 0 + F*F*3*r_SAR(1);
        BP_y = 0;
        BP_z = -F*F*r_SAR(3);
        a_3BP = [BP_x;BP_y;BP_z];
        a_3BP = SAR_N' * a_3BP;

        dU = dU + a_3BP;
    end

    % differential equations
       dx = [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3)];

end