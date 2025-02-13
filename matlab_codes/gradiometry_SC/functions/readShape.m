function [Nobj, X, Y, Z, R, D] = readShape(path)
%%                    READ SHAPE FUCNTIONS                               %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 12/27/2023                                                      %
%                                                                         %
%   Description: Read shape specified in the path directory.              %
%                                                                         %
%   Input: path: object shape directory                                   %
%   Output: Nobj: number of spheres                                       %
%           Xvec: X coordinates for object i                              %
%           Yvec: Y coordinates for object i                              %
%           Zvec: Z coordinates for object i                              %
%           Rvec: R size for object i                                     %
%           Dvec: Density for object i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define matrices
    Nmax = 5000;
    X = ones(Nmax, 1) * NaN;
    Y = ones(Nmax, 1) * NaN;
    Z = ones(Nmax, 1) * NaN;
    R = ones(Nmax, 1) * NaN;
    D = ones(Nmax, 1) * NaN;

    % scan file
    fid = fopen(path,'rt');
    j = 1;
    while true
      thisline = fgetl(fid);
      if ~ischar(thisline); break; end  %end of file
      ind1 = strfind(thisline, '--');
      ind2 = strfind(thisline, 'X');
      ind3 = strfind(thisline, 'SC');

      if(isempty(ind1) && isempty(ind2) && isempty(ind3))
          val = str2num(thisline);
          Nobj = val(1);
          X(j, 1) = val(2);
          Y(j, 1) = val(3);
          Z(j, 1) = val(4);
          R(j, 1) = val(5);
          D(j, 1) = val(6);
         
          j =  j + 1;
      end
    end
      fclose(fid);
      nanVal = isnan(X);
      ind = find(nanVal == 0);
      X = X(ind);
      Y = Y(ind);
      Z = Z(ind);
      R = R(ind);
      D = D(ind);
end

