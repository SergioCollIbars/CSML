function [ACI_ECEF] = read_ECEF2ITRF_mat(R)
    Nt = length(R(:, 1));
   
    % output matrix
    ACI_ECEF  = ones(3*Nt, 3) * NaN;
    for j = 1:Nt
        b1 = R(j, 2); b2 = R(j, 3); b3 = R(j, 4); b0 = R(j, 5);
       
        % construct matrix
        mat11 = b0^2 + b1^2 - b2^2 - b3^2;
        mat12 = 2*(b1*b2 + b0*b3);
        mat13 = 2*(b1*b3 - b0*b2);
        mat21 = 2*(b1*b2 - b0*b3);
        mat22 = b0^2 - b1^2 + b2^2 - b3^2;
        mat23 = 2*(b2*b3 + b0*b1);
        mat31 = 2*(b1*b3 + b0*b2);
        mat32 = 2*(b2*b3 - b0*b1);
        mat33 = b0^2 - b1^2 - b2^2 + b3^2;

        % ensamble matrix
        maxPos = 3 * j; minPos = maxPos - 2;
        ACI_ECEF(minPos:maxPos, :) = [mat11, mat12, mat13;mat21, mat22, mat23;...
            mat31, mat32, mat33];
    end
end

