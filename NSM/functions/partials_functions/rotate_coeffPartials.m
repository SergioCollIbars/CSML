function [Hc_BODY] = rotate_coeffPartials(Hc_ACI,BN)
    % Description: Rotate a (9 x Nc+Ns) matrix from the ACI frame to the 
    % body frame, BN rotation matrix.
    % Author: Sergio Coll Ibars
    % Date: 07/03/2025
    % Input:    Hc_ACI = coefficient partials matrix in ACI frame.
    %           BN = ACI to BODY frame rotation.
    % Output:   Hc_BODY = coefficient partials matrix in BODY frame.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define output
    Hc_BODY = Hc_ACI.* NaN;
    h = Hc_ACI;

    % pre-compute the rotation transpose
    BNt = BN';
    % fill-out rest of the terms
    for j = 1:length(Hc_ACI(1, :))
        TN = [h(1, j), h(4, j), h(7, j);h(2, j),h(5, j),h(8, j);...
                  h(3, j),h(6, j), h(9, j)];
        TB = BN * TN * BNt; % rotate to body frame
        Hc_BODY(:, j) = [TB(1,1);TB(2,1);TB(3,1);TB(1,2);TB(2,2);TB(3,2);...
        TB(1,3);TB(2,3);TB(3,3)];
    end
end

