function [He, b, S] = H_extraParam(varObj, k, At)
    %%                  H MATRIX EXTRA PARAMERTERS FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 11/05/2023
    %
    %   Description: This function returns the visibility matrix and bias
    %   vector for different extra parameters modes.
    %
    %   Input:
    %       varObj: estimation class object
    %       k : iteration number
    %       At: time increment
    %  
    %   Output: 
    %        H: visibility matrix
    %        b: bias matrix
    %        S: process noise covariance matrix
    % --------------------------------------------------------------------%


    if(varObj.extraParam == 1)
         % bias 
         Nc = varObj.Nc;
         Ns = varObj.Ns;
         Np = varObj.Np;

         init = Nc + Ns + 1;
         final = init + 11;
         val =  varObj.Xfilter(init:final, k);
         b = [val(1), val(3), val(5);...
             0, val(7), val(9);...
             0, 0, val(11)];
         He = [-1, zeros(1, 11);...
             zeros(1, 2), -1, zeros(1, 9);...
             zeros(1, 4), -1, zeros(1, 7);...
             zeros(1, 6), -1, zeros(1, 5);...
             zeros(1, 8), -1, zeros(1, 3);...
             zeros(1, 10), -1, 0];

         % process noise
         q = 1E-25;
         SE_left = zeros(12, Nc + Ns);
         SE_up = SE_left';
         At = 10.0012;
         se =  q*[At^3/3, At^2/2;At^2/2, At];
         se = [q , 0; 0, 0];
         SE = blkdiag(se, se, se, se, se, se);
         S = [zeros(Nc+Ns, Nc+Ns), SE_up;...
            SE_left, SE];
         S = zeros(Np, Np);
    elseif(varObj.extraParam == 2)
        S = 0;
        b = 0;
        
        Sb = 1E-10;
        Sd = 1E-10;
        He = [-1*Sb, -k*At*Sd, zeros(1, 10);...
             zeros(1, 2), -1*Sb, -k*At*Sd, zeros(1, 8);...
             zeros(1, 4), -1*Sb, -k*At*Sd, zeros(1, 6);...
             zeros(1, 6), -1*Sb, -k*At*Sd, zeros(1, 4);...
             zeros(1, 8), -1*Sb, -k*At*Sd, zeros(1, 2);...
             zeros(1, 10), -1*Sb, -k*At*Sd];
    else
        % no extra parameters
        b = 0;
        He = 0;
        S = zeros(varObj.Np, varObj.Np);
    end

end

