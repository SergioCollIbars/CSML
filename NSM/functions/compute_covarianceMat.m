function [P] = compute_covarianceMat(A)
    % Description: compute the covariance matrix given the information
    % matrix A.

    % Cholesky decomposiiton
    [R, p] = chol(A);
    if(p == 0) % matrix is positive definite
        P = inv(R) * inv(R');
    else        % not positive definite. Use pseudo-inverse
        P = pinv(A);
    end

    % check  covariance inversion quality
    I = eye(size(A));
    residual = norm(A * P - I, 'fro');  % Frobenius norm of residual
    disp([' Residual inversion = ', num2str(residual)]);
end
