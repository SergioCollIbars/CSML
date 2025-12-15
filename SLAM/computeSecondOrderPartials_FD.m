function [H2_mat] = computeSecondOrderPartials_FD(n_max, normalized, Cnm, Snm, Re, GM, ...
                r0_ACI, ACAF_ACI, ACAF_B)
%--------------------------------------------------------------------------
% Computes second-order partial derivatives (Hessian tensor) of measurements
% using finite differences of the FIRST-order Jacobian.
%
% INPUTS:
%   jacobianFun : function handle -> H = jacobianFun(r)
%                 must return 6x3 matrix
%   r0          : 3x1 nominal position
%   h           : finite difference step (e.g. 1e-3 to 1e-1 meters)
%
% OUTPUT:
%   H2          : 6 x 3 x 3 tensor
%                 H2(:,j,m) = d^2 y / (dr_j dr_m)
%--------------------------------------------------------------------------

    nMeas  = 6;
    nState = 3;

    H2 = zeros(nMeas, nState, nState);

    h = 1E-2;
    
    r0 = r0_ACI;

    % Loop over second derivative directions
    for m = 1:nState

        dr = zeros(3,1);
        dr(m) = h;

        % First-order Jacobians at perturbed points
        [Hpos_p] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, ...
            GM, r0 + dr, ACAF_ACI, ACAF_B);
        [Hpos_m] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, ...
            GM, r0 - dr, ACAF_ACI, ACAF_B);

        Hp       = [Hpos_p(1:3, :); Hpos_p(5:6, :);Hpos_p(9, :)];
        Hm       = [Hpos_m(1:3, :); Hpos_m(5:6, :);Hpos_m(9, :)];

        % Second-order derivative via central difference
        % ∂/∂r_m ( ∂y/∂r ) ≈ (Hp - Hm) / (2h)
        dHdm = (Hp - Hm)./(2*h);    % 6x3

        % Store each column as the second derivative wrt (j,m)
        for j = 1:nState
            H2(:, j, m) = dHdm(:, j);
        end
    end

    % order Hessian
    Hqq = compressHessianToQuadraticForm(H2);
    H2_mat = 0.5 * Hqq;
end
