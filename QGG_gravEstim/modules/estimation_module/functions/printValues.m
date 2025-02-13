function printValues(X, n_max, Nxc, Nxs)
   %%                    PRINT VALUES FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 25/03/2023
    %
    %   Description: This function print in the terminal the estimation
    %   values in a nice way.
    %
    %   Input:
    %       X: estimation vector
    %       n_max: max harmonics estimated order
    %
    %   Output: 
    %       N/A
    %
    % --------------------------------------------------------------------%
    format long g;

    disp("State vector: ");
    disp("C_11 = " + string(X(1)));
% %     for j = 2:n_max
% %         for i = 1:j+1
% %             disp("C_" + string(j) + string(i-1) + " = " + string(X(j + i - 1)));
% %         end
% %     end
% % 
% %     disp(" ");
% %     n = 2;
% %     i = 1;
% %     for j = Nxc+1:Nxc+Nxs
% %         val = X(j);
% %         disp("S_" + string(n) + string(i) + " = " + string(val));
% %         if(i < n)
% %             i = i + 1;
% %         else
% %             n = n + 1;
% %             i = 1;
% %         end
% %     end

    Xc = X(1:Nxc);
    Xs = X(1+Nxc:Nxc+Nxs);

    for n = 2:n_max
        c = -1;
        r = 0;
        for i  = 1:n
            c = c + i;
        end
        for j = c:c+n
            disp("C_" + string(n) + string(r) + " = " + string(Xc(j)));
            r = r + 1;
        end
    end
    for n = 2:n_max
        r = 0;
        c = 0;
        for i = 1:n-1
            c = c + i;
        end
        for j = c:c+n-1
            disp("S_" + string(n) + string(r) + " = " + string(Xs(j)));
            r = r + 1;
        end

    end

end

