function [U, Udot, Uddot] = grav_sphere(Nobj, Dvec, Xvec, Yvec, Zvec, ...
    Rvec, gradiometerPos)
%%                GRAVITATIONAL SPHERE FUCNTIONS                         %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 12/30/2023                                                      %
%                                                                         %
%   Description: Read shape specified in the path directory.              %
%                                                                         %                               
%   Input: Nobj: number of spheres                                        %
%           Dvec: j sphere density                                        %
%           Xvec: X coordinates for object i                              %
%           Yvec: Y coordinates for object i                              %
%           Zvec: Z coordinates for object i                              %
%           Rvec: R size for object i                                     %
%                                                                         %
%   Output:  U: gravity potential                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Gravirty contants
    G = 6.37e-11;
    
    % position vector
    r = gradiometerPos;

    % gravity potential, acceleration and 2nd order gravity tensor
    U = 0;
    Udot = zeros(3, 1);
    Uddot = zeros(3, 3);

    for j = 1:Nobj
        % Total mass
        M = Dvec(j) * 4/3 * pi * Rvec(j)^3;

        % Position
        rr = [Xvec(j); Yvec(j); Zvec(j)];
        s = vecnorm(r - rr);
        sv = r - rr;
        s3 = s^3;
        s5 = s^5;
        
        if(s >= Rvec(j)) % outside sphere j
            U = U - G * M / s;
            Udot = Udot + G*M/(s*s) * sv;
            Uddot = Uddot -G*M*[1/s3-3*sv(1)*sv(1)/s5,-3*sv(2)*sv(1)/s5, ...
                -3*sv(1)*sv(3)/s5;-3*sv(2)*sv(1)/s5,1/s3-3*sv(2)*sv(2)/s5, ...
                -3*sv(2)*sv(3)/s5;-3*sv(3)*sv(1)/s5,-3*sv(2)*sv(3)/s5, ...
                1/s3-3*sv(3)*sv(3)/s5];
        else % inside sphere j
            U = U + G * M / (2 * Rvec(j)^3) * (s*s - 3*Rvec(j)*Rvec(j));
            Udot = Udot + G*M/(Rvec(j)^3) * sv;
            Uddot = Uddot + G * M / (Rvec(j)^3) * eye(3,3);
        end
    end
end

