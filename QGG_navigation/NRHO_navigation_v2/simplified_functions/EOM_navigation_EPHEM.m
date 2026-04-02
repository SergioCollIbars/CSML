function [dx] = EOM_navigation_EPHEM(t, x, planetParams, C_mat, S_mat)
    %%                       EOM FUNCTION EPHEM
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 10/14/2025
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
    Ns = 12;

    % S/C parameters
    m = planetParams(11);   % [Kg]
    A = planetParams(12);   % [-]
    eta = 1.3;              % [-]

    % extract planet parameters (non-dimensional units)
    [GM2] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    [Re2] = cspice_bodvrd('MOON', 'RADII', 3); % get Sun Radius [Km]
    GM2 = GM2*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re2 = Re2*1E3./planetParams(2);

    [GM1] = cspice_bodvrd('EARTH', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
    [Re1] = cspice_bodvrd('EARTH', 'RADII', 3); % get Jupiter Radius [Km]
    GM1 = GM1*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re1 = Re1*1E3./planetParams(2);
    
    [GM3] = cspice_bodvrd('SUN', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    [Re3] = cspice_bodvrd('SUN', 'RADII', 3); % get Sun Radius [Km]
    GM3 = GM3*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re3 = Re3*1E3./planetParams(2);

    [GM4] = cspice_bodvrd('5', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
    [Re4] = cspice_bodvrd('JUPITER', 'RADII', 3); % get Jupiter Radius [Km]
    GM4 = GM4*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re4 = Re4*1E3./planetParams(2);

    % gravity computation params
    n_max      = planetParams(6);
    normalized = planetParams(7);
    
    % compute Earth position. ref: J2000
    target = 'EARTH';
    et = t./planetParams(3);      % Convert UTC time to ephemeris time
    et1 = t./planetParams(3) - 1; % Convert UTC time to ephemeris time
    et2 = t./planetParams(3) + 1; % Convert UTC time to ephemeris time
    ref = 'J2000';
    abcorr = 'NONE';
    observer = '3';  % Set the observer to the Earth-Moon barycenter

    [Estate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Epos  = Estate(1:3)*1E3/planetParams(2);                                           % [m]
    r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

    % compute Moon position. ref: J2000
    target = 'MOON';
    [Mstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Mpos  = Mstate(1:3)*1E3/planetParams(2);                                          % [m]
    r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

    % compute Sun position. ref: J2000
    target = 'SUN';
    [Sstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Spos = Sstate(1:3)*1E3/planetParams(2);                                            % [m]
    r3 = [x(1)-Spos(1);x(2)-Spos(2);x(3)-Spos(3)];                      % SC-Sun
    
    
    % EM barycenter acceleration. ref: J2000
    target = '3';
    [EMstate1, ~] = cspice_spkezr(target, et1, ref, abcorr, '0');       % [Km & Km/s]
    [EMstate2, ~] = cspice_spkezr(target, et2, ref, abcorr, '0');       % [Km & Km/s]

    Svel2 = EMstate2(4:6)*1E3;   % [m/s]
    Svel1 = EMstate1(4:6)*1E3;   % [m/s]
    
    
    At = (et2-et1);
    Acc_EM = (Svel2- Svel1)./At; % [m/s^2]
    Acc_EM = Acc_EM./(planetParams(2) * planetParams(3)^2); % [-]

    % compute Jupiter barycenter position. ref: J2000
    target = '5';
    [Jstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Jpos = Jstate(1:3)*1E3/planetParams(2);                             % [-]
    
    r4 = [x(1)-Jpos(1);x(2)-Jpos(2);x(3)-Jpos(3)];                      % SC-Sun

    % compute orientation
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et);

    % compute gravity acceleration
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, dU1, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                J2000_EARTH'*r1, Re1(1), GM1, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, dU2, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                J2000_MOON'*r2, Re2(1), GM2, ...
                                                normalized);
    Cmat3 = Cmat2;
    Smat3 = Smat2;
    Cmat3(2:end, :) = 0;
    Smat3(2:end, :) = 0;
    [~, dU3, ddU3] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                r3, Re3(1), GM3, ...
                                                normalized);

    [~, dU4, ddU4] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                r4, Re4(1), GM4, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    dU1  = J2000_EARTH  * dU1;
    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';

    dU2  = J2000_MOON  * dU2;
    ddU2 = J2000_MOON  * ddU2  * J2000_MOON';

    % compute SRP acceleration
    [aSRP, daSRP_dr, ~] = SRP(r3, eta, m, A);

    % total acceleration
    dU = dU1 + dU2 + dU3 + dU4 - Acc_EM + aSRP;

    % compute gravity position partials
    T = ddU1 + ddU2 + ddU3 + ddU4;
    
    % compute Jacobian
    J = compute_jacobian(T, daSRP_dr);

    % STM value 
    PHI = reshape(x(Ns+1:Ns + Ns*Ns), [Ns, Ns]);

    PHI_dot = J * PHI;  

    % bias acceleration
    b_dot = zeros(6, 1);

    % differential equations
   dx =  [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3);
          b_dot; reshape(PHI_dot, [Ns*Ns, 1])];
end

