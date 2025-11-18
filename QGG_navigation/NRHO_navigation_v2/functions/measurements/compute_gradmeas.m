function [T] = compute_gradmeas(t, x, planetParams, poleParams,...
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

    % get body data
    GM = planetParams(1);
    Re = planetParams(2);
    n_max = planetParams(3);
    Normalized = planetParams(4);
    
    % get pole parameters data
    W = poleParams(1);
    W0 = poleParams(2);
    RA = poleParams(3);
    DEC = poleParams(4);

    % gradiometer measurements
    T = zeros(9, length(t));

    for j = 1:length(t)
        if(enviroment == "2BP" || enviroment == "F2BP")
            % ACI position vector
            r_ACI = [x(j, 1); x(j, 2); x(j, 3)];
        
            % Inertial to ACAF rotation 
            Wt = W * t(j) + W0;
            ACAF_N = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        
            % ACAF position vector
            r_ACAF = ACAF_N * r_ACI;
            
            % compute S/C acc
            [~, ~, ddU] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                        r_ACAF, Re, GM, ...
                                                        Normalized);
            ddU = ACAF_N' * ddU * ACAF_N;

            % store measurements
            T(:, j) = reshape(ddU, [9, 1]);
    
        elseif(enviroment == "CR3BP" || enviroment == "3BP") % non-dimensional units
                mu = planetParams(1);
                
                r1 = sqrt((x(j, 1) + mu)^2 + x(j, 2)^2 + x(j, 3)^2);
                r2 = sqrt((x(j, 1) + mu - 1)^2 + x(j, 2)^2 + x(j, 3)^2);
   
    
                % compute SOGT
                ddU = gradmeas_rotFrame(mu, x(j, 1), x(j, 2), x(j, 3), ...
                    r1, r2);

                 % store measurements
                T(:, j) = reshape(ddU, [9, 1]);
        end
    end
end
