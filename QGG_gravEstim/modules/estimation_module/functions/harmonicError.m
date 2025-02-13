function [CS_err] = harmonicError(X, C_mat, S_mat, GM, Nxc, Nxs)
   %%                    PRINT VALUES FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 25/03/2023
    %
    %   Description: This function computes the harmonics error.
    %
    %   Input:
    %       X: estimation vector
    %       C_mat: C coefficient harmonics matrix
    %       S_mat: S coefficient harmonics matrix
    %
    %   Output: 
    %       CS_err: relative error matrix
    %
    % --------------------------------------------------------------------%
   
    % number of estimated pararms
    Nx = length(X);

    % error matrix definition
    CS_err = zeros(Nx, 1);

    CS_err(1) = abs(GM - X(1));
    
    % C_mat err
    n = 3;
    i = 1;
    for j = 2:Nxc
        CS_err(j) = abs(C_mat(n, i) - X(j));
        if(i < n)
            i = i + 1;
        else
            n = n + 1;
            i = 1;
        end
    end
    n = 3;
    i = 2;
    for j = Nxc+1:Nxc+Nxs
        CS_err(j) = abs(S_mat(n, i) - X(j));
        if(i < n)
            i = i + 1;
        else
            n = n + 1;
            i = 2;
        end
    end


end

