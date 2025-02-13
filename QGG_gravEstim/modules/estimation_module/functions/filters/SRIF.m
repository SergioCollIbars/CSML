function [R, b] = SRIF(dY, Hi, R, b, PHI_ji)
    %%                  SQUARE ROOT INFORMATION FILTER FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/01/2023
    %
    %   Description: Compute the normal equation in the information domain
    %                Allows to bigger numerical precission.
    %
    %   Input: 
    %       dY: prefit valies @ current time
    %       Hi: vidibility matrix @ current time
    %       R: UT information matrix @ t: i -1
    %       b: state information domain @ t: i-1
    %       PHI_ji: STM from t: i-1 to t: i

    %   Output:
    %       R: meas update UT information matrix @ t: i
    %       b: meas update state information domain @ t: i
    % --------------------------------------------------------------------%
  
    
    % values time update
    b_bar = b;
    R_bar = R * PHI_ji;
    
    % UT time update matrix
    n = length(R(1, :));
    A = [R_bar,b_bar];

    A = householderT(A);
    R_bar = A(1:n, 1:n);
    b_bar = A(1:n, n+1);

    if(isnan(dY))
        % update meas
        b = b_bar;
        R = R_bar;
    else
        % UT meas update
        A = [R_bar, b_bar; ...
             Hi, dY]; 

        A = householderT(A);
        
        % meas update values
        R = A(1:n, 1:n);
        b = A(1:n, n+1);
    end
end