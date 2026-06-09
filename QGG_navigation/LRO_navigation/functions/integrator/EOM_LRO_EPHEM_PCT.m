function [dx] = EOM_LRO_EPHEM_PCT(t, x, planetParams, C_mat, S_mat, n_maxM, n_maxE, tmin, tmax)
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
    
    %% Progress printing
    pct = 100 * (t - tmin)/(tmax - tmin);
    fprintf('\b\b\b\b%3.0f%%', pct);

    %% Constants

    Nx = 6;

    GM_E = planetParams(1);
    GM_M = planetParams(2);
    R_E  = planetParams(3);
    R_M  = planetParams(4);
    normalized = planetParams(6);

    GM_S = planetParams(7);
    R_S  = planetParams(9);

    M   = 220;
    A   = 2.8;
    eta = 2.5;

    r_sc = x(1:3);
    v_sc = x(4:6);

    %% SPICE states, Moon-centered J2000

    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';

    Estate = cspice_spkezr('EARTH', t, ref, abcorr, observer);
    Sstate = cspice_spkezr('SUN',   t, ref, abcorr, observer);

    Epos = Estate(1:3) * 1e3;   % Earth wrt Moon, J2000 [m]
    Spos = Sstate(1:3) * 1e3;   % Sun wrt Moon, J2000 [m]

    r1 = r_sc - Epos;           % spacecraft wrt Earth
    r2 = r_sc;                  % spacecraft wrt Moon
    r3 = r_sc - Spos;           % spacecraft wrt Sun

    r_ME = -Epos;               % Moon wrt Earth
    r_MS = -Spos;               % Moon wrt Sun

    %% Frame rotations

    J2000_EARTH = cspice_pxform('IAU_EARTH', 'J2000', t);
    J2000_MOON  = cspice_pxform('MOON_PA',   'J2000', t);

    %% SRP

    F = shadow_function(R_S, R_M, Spos, -r2);
    [aSRP, daSRP_dr, ~] = SRP(r3, eta, M, A);

    %% Earth gravity

    Cmat_E = C_mat{1};
    Smat_E = S_mat{1};

    r1_E = J2000_EARTH.' * r1;

    [~, dU1_E, ddU1_E] = potentialGradient_nm_mex( ...
        Cmat_E, Smat_E, n_maxE, r1_E, R_E, GM_E, normalized);

    dU1  = J2000_EARTH * dU1_E;
    ddU1 = J2000_EARTH * ddU1_E * J2000_EARTH.';

    % Point-mass indirect Earth acceleration
    a_tidal_E = GM_E * r_ME / norm(r_ME)^3;

    %% Moon gravity

    Cmat_M = C_mat{2};
    Smat_M = S_mat{2};

    r2_M = J2000_MOON.' * r2;

    [~, dU2_M, ddU2_M] = potentialGradient_nm_mex( ...
        Cmat_M, Smat_M, n_maxM, r2_M, R_M, GM_M, normalized);

    dU2  = J2000_MOON * dU2_M;
    ddU2 = J2000_MOON * ddU2_M * J2000_MOON.';

    %% Sun point-mass gravity and indirect acceleration

    [dU3, ddU3] = point_mass_accel_partials(r3, GM_S);

    a_tidal_S = GM_S * r_MS / norm(r_MS)^3;

    %% Lunar Albedo force
    [a_alb, da_alb_dr] = lunar_albedo_accel(r2, r_MS, 1.3, A, M, 0.12);

    %% Relativistic corrections
    [a_rel, da_rel_dr, da_rel_dv] = relativistic_accel(r2, v_sc, GM_M);

    %% Total acceleration and partials

    a_total = dU2 + dU1 + dU3 + a_tidal_E + a_tidal_S + F*aSRP + ...
        a_alb + a_rel;

    T = ddU2 + ddU1 + ddU3 + F*daSRP_dr + da_rel_dr + da_alb_dr;

    %% STM propagation

    PHI = reshape(x(Nx+1:Nx+Nx*Nx), [Nx, Nx]);

    J = zeros(Nx, Nx);
    J(1:3, 4:6) = eye(3);
    J(4:6, 1:3) = T;
    J(4:6, 4:6) = da_rel_dv;

    PHI_dot = J * PHI;

    %% Differential equations

    dx = [v_sc;
          a_total;
          PHI_dot(:)];
end