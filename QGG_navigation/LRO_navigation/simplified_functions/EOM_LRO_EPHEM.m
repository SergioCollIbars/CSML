function [dx] = EOM_LRO_EPHEM(t, x, planetParams, C_mat, S_mat)
    %%                       EOM FUNCTION EPHEM
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 09/24/2025
    %
    %   Description: Compute the EoM using the Ephemerides model in
    %   cislunar space.
    %
    %   Input:
    %       t: time vector
    %       x: state vector [r1, r2, r3, v1, v2, v3, eta, bias vector]'
    %       planetParams: planet parameters 
    %                     [GM, Re, nmax, normalized]
    %       poleParams: pole parameters
    %                   [W, W0, RA, DEC]
    %       Cmat: SH C coefficients
    %       Smat: SM S coefficients
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%
    
    % number of states
    Ns = 6;

    % extract planet parameters (non-dimensional units)
    [GM2] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Moon [km^3/s^2]
% %     [Re2] = cspice_bodvrd('MOON', 'RADII', 3); % get Moon Radius [Km]
    GM2 = GM2*1E9;                             % [m^3 s^-2]
% %     Re2 = Re2*1E3;                             % [m]
    
    Re2 = 1738*1E3;
% %     GM2 = 4902.80007*1E9;

    [GM1] = cspice_bodvrd('EARTH', 'GM', 1);    % Get GM for the Earth [km^3/s^2]
    [Re1] = cspice_bodvrd('EARTH', 'RADII', 3); % get Earth Radius [Km]
    GM1 = GM1*1E9;
    Re1 = Re1*1E3;
    
% %     Re1 = 6378.1363*1E3;
% %     GM1 = 398600.4354*1E9;
    
    [GM3] = cspice_bodvrd('SUN', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    GM3 = GM3*1E9;

% %     GM3 = 132712440041.94 * 1E9;


    [GM4] = cspice_bodvrd('5', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
    GM4 = GM4*1E9;

    % gravity computation params
    normalized = planetParams(7);
    
    % compute Earth position. ref: J2000
    target = 'EARTH';
    et = t;                         % Convert UTC time to ephemeris time

    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';             % Set the observer to the MOON

    [Estate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Epos  = Estate(1:3)*1E3;                                            % [m]
    r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

    % compute Moon position. ref: MOON
    target = 'MOON';
    [Mstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Mpos  = Mstate(1:3)*1E3;                                            % [m]
    r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

    % compute Sun position. ref: MOON
    target = 'SUN';
    [Sstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Spos = Sstate(1:3)*1E3;                                             % [m]
    r3 = [x(1)-Spos(1);x(2)-Spos(2);x(3)-Spos(3)];                      % SC-Sun

     target = '5';
    [Jstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Jpos = Jstate(1:3)*1E3;
    
    r4 = [x(1)-Jpos(1);x(2)-Jpos(2);x(3)-Jpos(3)];                      % SC-Sun

    % compute orientation
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et);

    % compute gravity acceleration
    r_ME = Mpos - Epos;
    r_MS = Mpos - Spos;
    r_MJ = Mpos - Jpos;

    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, dU1, ddU1] = potentialGradient_nm(Cmat1, Smat1, 0, ...
                                                J2000_EARTH'*r1, Re1(1), GM1, ...
                                                normalized);

    [~, dU1_T, ~] = potentialGradient_nm(Cmat1, Smat1, 0, ...
                                                J2000_EARTH'*r_ME, Re1(1), GM1, ...
                                                normalized);

    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, dU2, ddU2] = potentialGradient_nm(Cmat2, Smat2, 100, ...
                                                J2000_MOON'*r2, Re2(1), GM2, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    dU1  = J2000_EARTH  * dU1;

    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';

    dU2  = J2000_MOON   * dU2;
    ddU2 = J2000_MOON   * ddU2  * J2000_MOON';

    dU3  = GM3 * (r3 / (vecnorm(r3)^3));     % Sun acceleration. Point mass

    dU4  = GM4 * (r4 / (vecnorm(r4)^3));     % Jupiter acceleration. Point mass

    % tidial acceleration
    a_tidial_E = - J2000_EARTH  * dU1_T;
    a_tidial_S = - GM3 * r_MS / (vecnorm(r_MS)^3);
    a_tidial_J = - GM4 * r_MJ / (vecnorm(r_MJ)^3);
    
    % total acceleration
    dU = dU2 + dU1 +  dU3 + dU4 + a_tidial_E + a_tidial_S + a_tidial_J;

    % compute gravity position partials
    T = ddU2 + ddU1;
    
    % compute Jacobian
    J = compute_jacobian(T);

    % STM value 
    PHI = reshape(x(Ns+1:Ns + Ns*Ns), [Ns, Ns]);

    PHI_dot = J * PHI;  

    % differential equations
   dx =  [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3);
          reshape(PHI_dot, [Ns*Ns, 1])];
end
