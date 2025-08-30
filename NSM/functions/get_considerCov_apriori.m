function [Mxx_bar, Mxc_bar, Mcc_bar] = get_considerCov_apriori(P0, Pc, Pxc)
    %%                    GET APRIORI CONSIDER COV VALUES
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 28/06/2024
    %
    %   Description: This function computes the apriori values for the
    %   consider covariance
    %
    %   Input:
    %       P0:
    %
    %   Output:
    %       Mxx_bar:  
    % --------------------------------------------------------------------%

    Mcc_bar = pinv(Pc - Pxc' * pinv(P0) * Pxc);
    Mxx_bar = pinv(P0 - Pxc * pinv(Pc) * Pxc');
    Mxc_bar = -Mxx_bar * Pxc * pinv(Pc);
end

