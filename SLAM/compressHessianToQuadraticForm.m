function Hqq = compressHessianToQuadraticForm(H2)
%--------------------------------------------------------------------------
% Converts full Hessian tensor (6x3x3) into a 6x6 quadratic-form matrix
% Columns correspond to:
% [x^2, x*y, x*z, y^2, y*z, z^2]
%
% INPUT:
%   H2  : 6 x 3 x 3  tensor, H2(i,j,m) = d²y_i/(dr_j dr_m)
%
% OUTPUT:
%   Hqq : 6 x 6 matrix
%         epsHO = 0.5 * Hqq * [dx^2; dx*dy; dx*dz; dy^2; dy*dz; dz^2]
%--------------------------------------------------------------------------

    Hqq = zeros(6,6);

    % x^2
    Hqq(:,1) = H2(:,1,1);

    % x*y (sum of mixed partials)
    Hqq(:,2) = H2(:,1,2) + H2(:,2,1);

    % x*z
    Hqq(:,3) = H2(:,1,3) + H2(:,3,1);

    % y^2
    Hqq(:,4) = H2(:,2,2);

    % y*z
    Hqq(:,5) = H2(:,2,3) + H2(:,3,2);

    % z^2
    Hqq(:,6) = H2(:,3,3);
end
