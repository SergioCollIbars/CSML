function [ACAF_ACI] = compute_planetAtt(poleParams , t, saveData, ACI_ECEF)
    Nt = length(t);
    ACAF_ACI = ones(3*Nt, 3) * NaN;
    if(saveData == 1)
        for j = 1:Nt
            maxPos = 3 * j; minPos = maxPos - 2;
            ACAF_ACI(minPos:maxPos, :) = ACI_ECEF(minPos:maxPos, :)';
        end
    else
        W = poleParams(1); W0 = poleParams(2); RA = poleParams(3); 
        DEC = poleParams(4);
        for j =  1:Nt
            Wt = W0 + W * t(j);
            R =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            maxPos = 3*j; minPos = maxPos - 2;
            ACAF_ACI(minPos:maxPos, :) = R;
        end
    end
end