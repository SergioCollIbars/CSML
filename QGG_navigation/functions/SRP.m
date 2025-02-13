function [aSRP, daSRP_dr, daSRP_dEta] = SRP(rs, eta, m, A, planetParams)
    %%                   SRP ACCELERATION FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 01/20/2023
    %
    %   Description: This function computes the SRP acceleration at certain
    %   time.
    %
    %   Input:
    %       rs: SC position to sun inerital frame. ECI frame
    %       eta: scale factor
    %       A: S/C reference area [-]
    %       m: S/C reference mass [kg]
    %
    %   Output: 
    %       aSRP: SRP acceleration. ACI frame
    %       daSRP_dr: partial respect to inertial position
    %       daSRP_dEta: partial respect to eta
    % --------------------------------------------------------------------%

    % SRP parameters
    Gamma = 5.67E-8 / (planetParams(3)^3);                % [Kg/K^4]
    c = 2.99792458E8 / (planetParams(2)*planetParams(3)); % [-] 
    Rs = 6.96E8 / planetParams(2);                        % [-]
    Ts = 5778;                                            % [K]
    Cr = 1;

    % SRP acc
    rsn = vecnorm(rs);
    rsu  = rs./rsn;
    P = (Gamma * Rs^2 * Ts^4 )/ (c * m);

    aSRP = eta * ( P * Cr * A) * rs / rsn^3;

    % SPR partial respect to inertial position
    daSRP_dr = eta * Gamma * Ts^4 * Rs^2 * Cr * A / (m * c * rsn^3) * ...
        (eye(3) - 3 * (rs * rs') / (rsn^2));

    % SRP partial respect to eta
    daSRP_dEta = ( P * Cr * A) * rs / rsn^3;
end

