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

    Mcc_bar = inv(Pc - Pxc' * inv(P0) * Pxc);
    Mxx_bar = inv(P0 - Pxc * inv(Pc) * Pxc');
    Mxc_bar = -Mxx_bar * Pxc * inv(Pc);
end

