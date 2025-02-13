function [dx] = EOM_navigation(t, x, planetParams, poleParams,...
    C_mat, S_mat, enviroment, consider_cov, process_noise, augmented_st)
    %%                          EoM FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 31/10/2022
    %
    %   Description: This function defines the equation of motion (EoM) for
    %   the orbital problem
    %
    %   Input:
    %       t: time vector
    %       x: state vector [r1, r2, r3, v1, v2, v3]'
    %       planetParams: planet parameters 
    %                     [GM, Re, nmax, normalized]
    %       poleParams: pole parameters
    %                   [W, W0, RA, DEC]
    %       Cmat: SH C coefficients
    %       Smat: SM S coefficients
    %       enviroment: 2BP dynamics or CR3BP
    %
    %   Output:
    %       dx:  diferential equation matrix
    % --------------------------------------------------------------------%

    % S/C parameters
    AU = 1.496e+11 / planetParams(2);   % [-]
    eta = planetParams(10);
    m = planetParams(11);
    A = planetParams(12);

    if(enviroment == "CR3BP" || enviroment == "FCR3BP")  % non-dimensional units
        % compute motion of the primaries. J2000 frame
        M = t;
        mu = planetParams(1);
            
        r1 = [x(1)+mu*cos(M);x(2)+mu*sin(M);x(3)];              % SC-Earth
        r2 = [x(1)+(mu - 1)*cos(M); x(2)+(mu-1)*sin(M); x(3)];  % SC-Moon
        r3 = [x(1)-AU*cos(M+pi/4); x(2)-AU*sin(M+pi/4); x(3)];  % SC-Sun

        % rotation from Earth-Moon planet to J2000
        i_EM = deg2rad(0);   % [rad]
        EM_N = rotationMatrix(0, 0, i_EM, [1, 1, 1]);

        % extract planet parameters (non-dimensional units)
        GM1 = planetParams(8)./(planetParams(2)^3 * planetParams(3)^2);
        GM2 = planetParams(9)./(planetParams(2)^3 * planetParams(3)^2);
        Re1 = planetParams(4)./planetParams(2);
        Re2 = planetParams(5)./planetParams(2);

        n_max      = planetParams(6);
        normalized = planetParams(7);

        % extract rotation parameters
        RA_E = poleParams(1);            % RA Earth [rad]
        DEC_E = poleParams(2);           % DEC Earth [rad]
        W0_E = poleParams(3);            % prime meridian Earth [rad]
        WDot_E = poleParams(4);          % ang. velocity Earth [rad/s]

        RA_M = poleParams(5);            % RA Moon [rad]
        DEC_M = poleParams(6);           % DEC Moon [rad]
        W0_M = poleParams(7);            % prime meridian Moon [rad]
        WDot_M = poleParams(8);          % ang. velocity Moon [rad/s]

        Wt_E = WDot_E * t / planetParams(3) + W0_E;
        Wt_M = WDot_M * t / planetParams(3) + W0_M;
        ACAF1_N = rotationMatrix(pi/2 + RA_E, pi/2 - DEC_E, Wt_E, [3, 1, 3]);
        ACAF2_N = rotationMatrix(pi/2 + RA_M, pi/2 - DEC_M, Wt_M, [3, 1, 3]);
        
        ACAF1_EM = ACAF1_N * EM_N';
        ACAF2_EM = ACAF2_N * EM_N';

        % compute gravity acceleration
        Cmat1 = C_mat{1};
        Smat1 = S_mat{1};
        [~, dU1, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                    ACAF1_EM*r1, Re1, GM1, ...
                                                    normalized);
        Cmat2 = C_mat{2};
        Smat2 = S_mat{2};
        [~, dU2, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                    ACAF2_EM*r2, Re2, GM2, ...
                                                    normalized);
        % rotate back to inertial. Earth-Moon (EM) plane
        dU1  = ACAF1_EM' * dU1;
        dU2  = ACAF2_EM' * dU2;

        ddU1 = ACAF1_EM' * ddU1 * ACAF1_EM;
        ddU2 = ACAF2_EM' * ddU2 * ACAF2_EM;
        
        % compute SRP acceleration
        [aSRP, daSRP_dr, daSRP_deta] = SRP(r3, eta, m, A, planetParams);

        % total acceleration
        dU = dU1 + dU2 + aSRP;

        % compute SOGT
        T = ddU1 + ddU2;

    elseif(enviroment == "EPHEM" || enviroment == "EPHEMP") % use ephemerides model, non-dimensional model
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
        Epos  = Estate(1:3)*1E3/planetParams(2);
        r1    = [x(1)-Epos(1);x(2)-Epos(2);x(3)-Epos(3)];                   % SC-Earth

        % compute Moon position. ref: J2000
        target = 'MOON';
        [Mstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Mpos  = Mstate(1:3)*1E3/planetParams(2);
        r2    = [x(1)-Mpos(1);x(2)-Mpos(2);x(3)-Mpos(3)];                   % SC-Moon

        % compute Sun position. ref: J2000
        target = 'SUN';
        [Sstate, ~] = cspice_spkezr(target, et, ref, abcorr, observer);     % [Km & Km/s]
        Spos = Sstate(1:3)*1E3/planetParams(2);
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
        Jpos = Jstate(1:3)*1E3/planetParams(2);
        
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
        [aSRP, daSRP_dr, daSRP_deta] = SRP(r3, eta, m, A, planetParams);

        % total acceleration
        dU = dU1 + dU2 + dU3 + dU4 - Acc_EM + aSRP;

        % compute gravity position partials
        T = ddU1 + ddU2 + ddU3 + ddU4;

        ACAF1_EM = J2000_EARTH';
        ACAF2_EM = J2000_MOON';
    end

    if(process_noise{1})
        % state number
        Ns = 9;

        % get unmodeled acc
        w = x(7:9);
        Wdot = - process_noise{2} * w;

        % compute Jacobian
        A = [zeros(3, 3), eye(3,3); T, zeros(3, 3)];
        J = [A, [zeros(3,3);eye(3,3)];zeros(3,6),-process_noise{2}];
    else
        if(augmented_st)
            % states number include eta factor
            Ns = 7;
            
            % compute Jacobian
            J = [zeros(3, 3), eye(3,3), zeros(3, 1); ...
                T+daSRP_dr, zeros(3, 3), daSRP_deta; ...
                zeros(1,7)];
        else
            % states number
            Ns = 6;
            
            % compute Jacobian
            J = [zeros(3, 3), eye(3,3); T+daSRP_dr, zeros(3, 3)];
        end
    end
    
    % STM value 
    PHI = reshape(x(Ns+1:Ns + Ns*Ns), [Ns, Ns]);

    PHI_dot = J * PHI;

   % consider convariance. WARNING: needs to be modified to account for DMC
    if(consider_cov)
        Nc = length(x) - 42;
        SIGMA = reshape(x(43:end), [6, Nc/6]);
        [b1, ~] = potentialGradient_Cnm(n_max, ACAF1_EM*r1, Re1, GM1, ...
            ACAF1_EM', normalized);
        [b2, ~] = potentialGradient_Cnm(n_max, ACAF2_EM*r2, Re2, GM2, ...
            ACAF2_EM', normalized);
        B = [zeros(3, Nc/6); b1(:, 2:end), b2(:, 2:end)];
        SIGMA_dot = J * SIGMA + B;
    end


    % differential equations
   dx = [x(4);
      x(5);
      x(6);
      dU(1);
      dU(2);
      dU(3)];

   if(process_noise{1})
       dx(4:6) = dx(4:6) + w;
       dx = [dx; Wdot; reshape(PHI_dot, [Ns*Ns, 1])];
   else
       if(augmented_st)
        dx = [dx; 0; reshape(PHI_dot, [Ns*Ns, 1])];
       else
        dx = [dx; reshape(PHI_dot, [Ns*Ns, 1])];
       end
   end
   if(consider_cov)
       dx = [dx;reshape(SIGMA_dot, [Nc, 1])];
   end
end
