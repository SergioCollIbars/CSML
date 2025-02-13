function [SRP] = SRP(OrbitObj, rs)
    %%                          SRP FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/3/2023
    %
    %   Description: This function resurts the SRP acceleration for a given
    %   position vector
    %
    %   Input:
    %       rs: SC position vector respect to Sun. ECI frame
    %       OrbitObj: Orbit object
    %
    %   Output:
    %       SRP:  SRP acceleration. ECI frame
    % --------------------------------------------------------------------%
    rs = rs*1E-3;
    % Spacecraft SRP model values
    etaSRP = 1;
    Cr = 1.28;
    A_m = 1 / 100;
    A_m = A_m * 1E-6;
    
    % position vector
    rsn = vecnorm(rs);
    rsu = rs./rsn;

    % solar pressure
    %P = OrbitObj.gamma * OrbitObj.Rs^2 * OrbitObj.Ts^4 / (OrbitObj.c);
    P = 1E14;

    % SRP acceleration
    SRP = (etaSRP * P * Cr * A_m) / (rsn^2) * rsu;

    % convert to m
    SRP = SRP * 1E3;

end

