function [x] = trapzIntegrator(a, b, f, t)
    %%                TRAPEZOIDAL INTEGRATION FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 05/15/2023
    %
    %   Description: This function computes the integral value of a    
    %       function using the trapezoidal approxiamtion.
    %
    %   Input:
    %       a: under integral limit. Escalar value
    %       b: upper integral limit. Escalar value
    %       f: function to be integrated. Vector 3xNt size
    %       t: integration time vector. Vector 1xNt size
    %
    %   Output:
    %       x: integral value. Vector 3X1 size
    %
    % --------------------------------------------------------------------%

    x = (t(b) -t(a)) * (f(:, a) + f(:, b)) * 0.5;
end

