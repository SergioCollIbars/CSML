function [dx] = EOM_3BP(t, x, planetParams, poleParams,...
    C_mat, S_mat, enviroment)
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
    %       planetParams: planet parameters 
    %                     [GM, Re, nmax, normalized]
    %       poleParams: pole parameters
    %                   [W, W0, RA, DEC]
    %       Cmat: SH C coefficients
    %       Smat: SM S coefficients
    %       enviroment: 2BP dynamics or CR3BP
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%

    % Get orbit data
    GM = planetParams(1);
    Re = planetParams(2);
    n_max = planetParams(3);
    Normalized = planetParams(4);

    W = poleParams(1);
    W0 = poleParams(2);
    RA = poleParams(3);
    DEC = poleParams(4);

    if(enviroment == "2BP")

        % ACI position vector
        r_ACI = [x(1); x(2); x(3)];
    
        % Inertial to ACAF rotation 
        Wt = W * t + W0;
        ACAF_N = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % ACAF position vector
        r_ACAF = ACAF_N * r_ACI;
        
        % compute S/C acc
        [~, dU, ddU] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                    r_ACAF, Re, GM, ...
                                                    Normalized);
        % rotate from ECEF 2 ECI
        dU = ACAF_N' * dU;
        ddU = ACAF_N' * ddU * ACAF_N;
        
        J = [zeros(3, 3), eye(3, 3);ddU, zeros(3, 3)];

        PHI = reshape(x(7:end), [6, 6]);
        PHI_dot = J  * PHI;

    elseif(enviroment == "CR3BP" || enviroment == "3BP") % non-dimensional units
            dU = zeros(3, 1);
            
            mu = planetParams(1);
            
            r1 = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
            r2 = sqrt((x(1) + mu - 1)^2 + x(2)^2 + x(3)^2);
            
            dU(1) = 2*x(5) + x(1) - (1-mu)*(x(1)+mu)/(r1^3) - ...
                mu*(x(1)+mu-1)/(r2^3);
            dU(2) = -2*x(4) + x(2) - (1-mu)*x(2)/(r1^3) - mu*x(2)/(r2^3);
            dU(3) = -(1-mu)*x(3)/(r1^3) - mu*x(3)/(r2^3);

            PHI = reshape(x(7:end), [6, 6]);

            % compute SOGT
            T = gradmeas_rotFrame(mu, x(1), x(2), x(3),r1, r2);

            % compute Jacobian
            J = [zeros(3, 3), eye(3,3); T, [0,2,0;-2,0,0;0,0,0]];

            PHI_dot = J * PHI;
    elseif(enviroment == "CR3BP_inertial")
            mu = planetParams(1);   % [-]

            % SC acceleration
            res = x(1:3) - [-mu*cos(t);-mu*sin(t);0];  % [m]
            rms = x(1:3) - [(1-mu)*cos(t);(1-mu)*sin(t);0];  % [m]

            dU = -((1-mu)/(vecnorm(res)^3)*res + mu/(vecnorm(rms)^3)*rms);

            ddU = gradmeas_CR3BP_inertial(mu,x(1), x(2),...
                x(3), t);

            J = [zeros(3, 3), eye(3, 3); ddU, zeros(3, 3)];

            PHI = reshape(x(7:end), [6, 6]);
            PHI_dot = J  * PHI;
    end

    % differential equations
   dx = [x(4);
      x(5);
      x(6);
      dU(1);
      dU(2);
      dU(3);
      reshape(PHI_dot, [36, 1])];
end
