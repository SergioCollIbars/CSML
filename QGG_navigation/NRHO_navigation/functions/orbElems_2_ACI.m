function [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e)
    % compute initial state vectors (orbit frame)
    r0 = rho / (1 + e * cos(f)) * [cos(f);...
                                  sin(f);...
                                  0];
    v0 = sqrt(GM / rho) * [-sin(f);...
                          e + cos(f);...
                          0];
    
    % compute rotation matrix: Body to ACI
    [BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
    
    % rotate initial state vectors. {I, J, K} frame
    r0 = BN' * r0;
    v0 = BN' * v0;
end

