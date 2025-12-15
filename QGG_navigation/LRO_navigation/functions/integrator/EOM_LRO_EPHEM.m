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
    
    % number of states (pos + vel + bias)
    Nx = 6;
    
    % Planet constants
    GM_E = planetParams(1); GM_M = planetParams(2);
    R_E  = planetParams(3); R_M  = planetParams(4);
    
    % Get GM for the Sun [m^3/s^2]
    [GM_S] = cspice_bodvrd('SUN', 'GM', 1)*1E9;    

    % gravity computation params
    normalized = planetParams(6);
    n_max      = planetParams(5); 
    
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

    % compute orientation
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et);

    % compute gravity acceleration
    Cmat_E = C_mat{1};
    Smat_E = S_mat{1};
    [~, dU1, ddU1] = potentialGradient_nm(Cmat_E, Smat_E, n_max, ...
                                                J2000_EARTH'*r1, R_E, GM_E, ...
                                                normalized);
    Cmat_M = C_mat{2};
    Smat_M = S_mat{2};
    [~, dU2, ddU2] = potentialGradient_nm(Cmat_M, Smat_M, n_max, ...
                                                J2000_MOON'*r2, R_M, GM_M, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    dU1  = J2000_EARTH  * dU1;
    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';

    dU2  = J2000_MOON  * dU2;
    ddU2 = J2000_MOON  * ddU2  * J2000_MOON';

    % Tidial acceleration
    r_SS = r3;
    r_ME = Mpos - Epos;
    r_MS = Mpos - Spos;
    a_tidial_E = dU1 + GM_E * r_ME / (vecnorm(r_ME)^3);
    a_tidial_S = GM_S * (r_SS / (vecnorm(r_SS)^3) - r_MS / (vecnorm(r_MS)^3));
    
    % total acceleration
    dU = dU2 + a_tidial_E + a_tidial_S;

    % compute gravity position partials
    T = ddU2 + ddU1;
    
    % compute Jacobian
    J = compute_jacobian(T);

    % STM value 
    PHI = reshape(x(Nx+1:Nx+Nx*Nx), [Nx, Nx]);

    PHI_dot = J * PHI;  

    % differential equations
   dx =  [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3);
          reshape(PHI_dot, [Nx*Nx, 1])];
end
