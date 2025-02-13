function [a3BP] = acc3BP(OrbitObj, rsp, rcp)
    %%                          SRP FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/22/2023
    %
    %   Description: This function resurts the 3BP acceleration for a given
    %   position vector
    %
    %   Input:
    %       rsp: SC to Sun position. ECI frame
    %       rcp: Planet to Sun position. ECI frame
    %       OrbitObj: Orbit object
    %
    %   Output:
    %       3BP:  3BP acceleration. ECI frame
    % --------------------------------------------------------------------%
    % Sun GM
    GM = OrbitObj.GMs;

    % 3BP acceleration
    rspn = vecnorm(rsp);
    rcpn = vecnorm(rcp);

    a3BP = -GM * (rsp/(rspn^3) - rcp/(rcpn^3));
end

