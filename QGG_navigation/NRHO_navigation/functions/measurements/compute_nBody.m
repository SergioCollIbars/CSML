function [ddU] = compute_nBody(x, t, C_mat, S_mat, planetParams, posE, posM, posS)
    
    ref = 'J2000';
    abcorr = 'NONE';
    observer = '3';  % Set the observer to the Earth-Moon barycenter

    % extract planet parameters (non-dimensional units)
    GM1 = planetParams(8)./(planetParams(2)^3 * planetParams(3)^2);
    GM2 = planetParams(9)./(planetParams(2)^3 * planetParams(3)^2);
    Re1 = planetParams(4)./planetParams(2);
    Re2 = planetParams(5)./planetParams(2);
    
    [GM3] = cspice_bodvrd('SUN', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    [Re3] = cspice_bodvrd('SUN', 'RADII', 3); % get Sun Radius [Km]
    GM3 = GM3*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re3 = Re3*1E3./planetParams(2);

    [GM4] = cspice_bodvrd('5', 'GM', 1);    % Get GM for the Jupiter [km^3/s^2]
    [Re4] = cspice_bodvrd('JUPITER', 'RADII', 3); % get Jupiter Radius [Km]
    GM4 = GM4*1E9/(planetParams(2)^3 * planetParams(3)^2);
    Re4 = Re4*1E3./planetParams(2);

    % compute orientation
    et = t /planetParams(3);                            
    frame_to   = 'J2000';
    frame_from = 'IAU_EARTH';
    J2000_EARTH = cspice_pxform(frame_from, frame_to, et);

    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et);

    % gravity computation params
    n_max      = planetParams(6);
    normalized = planetParams(7);
    
    r1    = [x(1)-posE(1);x(2)-posE(2);x(3)-posE(3)];                   % SC-Earth
    r2    = [x(1)-posM(1);x(2)-posM(2);x(3)-posM(3)];                   % SC-Moon
    r3    = [x(1)-posS(1);x(2)-posS(2);x(3)-posS(3)];                   % SC-Sun

    % compute Jupiter barycenter position. ref: J2000
    target = '5';
    [Jstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
    Jpos = Jstate(1:3)*1E3/planetParams(2);
    
    r4 = [x(1)-Jpos(1);x(2)-Jpos(2);x(3)-Jpos(3)];                      % SC-Sun

    % compute gravity acceleration
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, ~, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                J2000_EARTH'*r1, Re1(1), GM1, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, ~, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                J2000_MOON'*r2, Re2(1), GM2, ...                                                
                                                normalized);

    Cmat3 = Cmat2;
    Smat3 = Smat2;
    Cmat3(2:end, :) = 0;
    Smat3(2:end, :) = 0;
    [~, ~, ddU3] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                r3, Re3(1), GM3, ...
                                                normalized);

    [~, ~, ddU4] = potentialGradient_nm(Cmat3, Smat3, 0, ...
                                                r4, Re4(1), GM4, ...
                                                normalized);

    % rotate back to inertial. Earth-Moon (EM) plane
    ddU1 = J2000_EARTH  * ddU1  * J2000_EARTH';
    ddU2 = J2000_MOON   * ddU2  * J2000_MOON';

    % compute gravity position partials
    ddU = ddU1 + ddU2 + ddU3 + ddU4;
end

