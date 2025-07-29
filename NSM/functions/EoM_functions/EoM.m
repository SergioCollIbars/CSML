function [dx] = EoM(t, x, Cnm_t, Snm_t, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, augmented)
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
    %       augmented: augmented to the C and S coefficients.
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

    if(augmented)
        [Nc, Ns, ~] = count_num_coeff(n_max);
        [Hacc, ~] = potentialGradient_Cnm(n_max, r_ACAF, Re, GM, ...
                        ACAF_ACI', normalized);
         Nc = Nc-1;
        J = [zeros(Nc, Nc), zeros(Nc, Ns), zeros(Nc, 3), zeros(Nc, 3);...
            zeros(Ns, Nc), zeros(Ns, Ns), zeros(Ns, 3), zeros(Ns, 3);...
            zeros(3, Nc), zeros(3, Ns), zeros(3, 3), eye(3,3);...
            Hacc(:, 2:end), ddU_ACI, zeros(3,3)];
         % acount for GM
% %         [Hacc, ~] = potentialGradient_Cnm(n_max, r_ACAF, Re, GM, ...
% %                         ACAF_ACI', normalized);

% %         [Hrot] = partials_acc_EulerAngles(ACAF_ACI, Cnm_t, Snm_t, n_max, Re, GM, r_ACI, normalized, ACAF_ACI);
% % 
% %         J = [zeros(Nc, Nc), zeros(Nc, Ns), zeros(Nc, 3), zeros(Nc, 6);...
% %             zeros(Ns, Nc), zeros(Ns, Ns), zeros(Ns, 3), zeros(Ns, 6);...
% %             zeros(3, Nc), zeros(3, Ns), zeros(3, 3), eye(3,3), zeros(3, 3);...
% %             Hacc(:, 1:end), ddU_ACI, zeros(3,3), Hrot;...
% %             zeros(3, Nc), zeros(3, Ns), zeros(3, 6), zeros(3,3)];

% %                 J = [zeros(Nc, Nc), zeros(Nc, Ns), zeros(Nc, 3), zeros(Nc, 3);...
% %             zeros(Ns, Nc), zeros(Ns, Ns), zeros(Ns, 3), zeros(Ns, 3);...
% %             zeros(3, Nc), zeros(3, Ns), zeros(3, 3), eye(3,3);...
% %             Hacc(:, 1:end), ddU_ACI, zeros(3,3)];
        Nx  = Nc + Ns + 6;
    else
        Nx  = 6;
    end

    PHI = reshape(x(7:end), [Nx, Nx]);
    PHI_dot = J  * PHI;

    % differential equations
       dx = [x(4);
          x(5);
          x(6);
          dU_ACI(1);
          dU_ACI(2);
          dU_ACI(3);
          reshape(PHI_dot, [Nx*Nx, 1])];
end

function [Nc, Ns, Ncs] = count_num_coeff(degree)
    % DESCRIPTION: count the number of zonal and sectorial / tesseral
    % coeffcients for a degree (n) value.
    Nc = 1;
    for k = 2:degree
        Nc = Nc + k + 1;
    end
    Ns = 0;
    for k = 2:degree
        Ns = Ns + k;
    end
    
    % total number of coefficients
    Ncs = Nc + Ns;
end


function [Hrot] = partials_acc_EulerAngles(ACAF_ACI, Cmat, Smat, n_max, Re, GM, r, normalized, ACAF_B)

    % output value
    Hrot = ones(3, 3) * NaN;
    ACI_ACAF = ACAF_ACI';

    eps = 1E-6;
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;

        Atpos = At./2;
        Atneg = - At./2; 

        [Rpos] = rotationMatrix(Atpos(1), Atpos(2), Atpos(3), [3, 2, 1]);
        [Rneg] = rotationMatrix(Atneg(1), Atneg(2), Atneg(3), [3, 2, 1]);

        [~, ddUpos, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*r, Re, GM, ...
                                                normalized);    % output in ACAF
        ddUpos = ACAF_B' * ddUpos; % S/C body frame
        ddUpos = Rpos * ddUpos;

        [~, ddUneg, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                ACI_ACAF'*r, Re, GM, ...
                                                normalized);    % output in ACAF
        ddUneg = ACAF_B' * ddUneg; % S/C body frame
        ddUneg = Rneg * ddUneg;

        H = (ddUpos - ddUneg)./(vecnorm(Atpos-Atneg));

        Hrot(:, j) = H;
    end
end
