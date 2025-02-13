function [Xhat, Y, fit] = solver_initCond_LS(t, v0, At, STM, gamma)
    % batch solver
    N0 = 2;
    Nf = length(t)-1;
    N = Nf - N0;
    Y = zeros(3, N);
    H = zeros(9, N);
    Ax = 0;
    Nx = 0;
    for j = N0:Nf
        PHI_next = reshape(STM(j+1, 1:end), [6,6]);
        PHI_prev = reshape(STM(j-1, 1:end), [6,6]);
        PHI_j    = reshape(STM(j, 1:end), [6,6]);
        alpha = PHI_next - PHI_prev;
        
        a11 = alpha(1:3, 1:3);
        a12 = alpha(1:3, 4:6);
        a21 = alpha(4:6, 1:3);
        a22 = alpha(4:6, 4:6);
    
        p11 = PHI_j(1:3, 1:3);
        p12 = PHI_j(1:3, 4:6);
        p21 = PHI_j(4:6, 1:3);
        p22 = PHI_j(4:6, 4:6);
    
        T = reshape(gamma(:, j), [3,3]);
    
        h = a22 - T*p12*(2*At);
        y = T*p11*v0*(2*At) - a21*v0;
        Y(:, j) = y;    % measurement vector
        H(:, j) = reshape(h, [9, 1]);    % measurement partial vector
    
        Ax = Ax + (h' * h);
        Nx = Nx + (h' * y);
    end
    Xhat = Ax\Nx;
    
    fit = zeros(3, N);
    for j = N0:Nf
        h = reshape(H(:, j), [3,3]);
        fit(:, j) = h*Xhat;
    end
end

