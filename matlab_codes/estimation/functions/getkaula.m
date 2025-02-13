function [Cnm0,Snm0] = getkaula(n_max, Nc, Ns, normalized)
    %%                   GET KAULA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 11/01/2023
    %
    %   Description: This function returns the kaula parameters for the
    %   specified SH degree.
    %
    %   Input:
    %       n_max: SH max degree
    %       Normalized: 1 for normalized coefficients, 0 otherwise
    %  
    %
    %   Output: 
    %        Cnm0: initial value for Kaula Cnm parameters
    %        Snm0: inital value for kaula Snm parameters
    % --------------------------------------------------------------------%


    % Kaula normalized bounds. Bennu
    KaulaN_z = 0.084;
    KaulaN_s = 0.026;
    alphaN_s = 2.01;
    alphaN_z = 2.08;
    
    % 1st row: zonals / 2nd row: sectorials and tesserals
    Kaula = zeros(2, n_max); 

    for n= 1:n_max
        if(normalized == 1)
            N_zonal = 1;
            N_sectorial = 1;
        else
            [N_zonal] = NormFactor(n, 0);
            [N_sectorial] = NormFactor(n, 1);
        end

        Kaula(1, n) = (KaulaN_z/n^alphaN_z) * N_zonal;
        Kaula(2, n) = (KaulaN_s/n^alphaN_s) * N_sectorial;
    end
    
    % create matrices
    Cnm0 = zeros(Nc, 1);
    Snm0 = zeros(Ns, 1);

    % asign coefficients
    Cnm0(1) = 1;
    
    m  = 0;
    n = 2;
    for j =2:Nc
        if(m == 0)
            Cnm0(j) = Kaula(1, n);
        else
            Cnm0(j) = Kaula(2, n);
        end
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
end

