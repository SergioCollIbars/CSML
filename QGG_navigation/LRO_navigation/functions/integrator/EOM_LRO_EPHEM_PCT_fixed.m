function [dx] = EOM_LRO_EPHEM_PCT_fixed(t, x, planetParams, C_mat, S_mat, n_maxM, n_maxE, tmin, tmax)
    %%                       EOM FUNCTION EPHEM
    % ------------------------------------------------------------------- %
    % Fixed version: Moon-centered J2000 dynamics with consistent vector
    % definitions, shadow geometry, albedo Sun vector, and third-body terms.
    %
    % State: x = [r_MSC; v_MSC; Phi(:)]
    %   r_MSC : Moon -> spacecraft, expressed in J2000 [m]
    %   v_MSC : spacecraft velocity wrt Moon, expressed in J2000 [m/s]
    % ------------------------------------------------------------------- %

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

    M   = 220;      % spacecraft mass [kg]
    A   = 2.8;      % effective area [m^2]
    % % eta = 1.75;      % SRP coefficient used by SRP()
    eta = 0;

    r_sc_M = x(1:3);    % Moon -> spacecraft, J2000 [m]
    v_sc_M = x(4:6);    % spacecraft wrt Moon, J2000 [m/s]

    %% SPICE states, Moon-centered J2000
    ref      = 'J2000';
    abcorr   = 'NONE';
    observer = 'MOON';

    Estate = cspice_spkezr('EARTH', t, ref, abcorr, observer);
    Sstate = cspice_spkezr('SUN',   t, ref, abcorr, observer);

    r_E_M = Estate(1:3) * 1e3;    % Moon -> Earth, J2000 [m]
    r_S_M = Sstate(1:3) * 1e3;    % Moon -> Sun,   J2000 [m]

    % Relative vectors with explicit origin convention
    r_SC_E = r_sc_M - r_E_M;      % Earth -> spacecraft [m]
    r_M_E  = -r_E_M;              % Earth -> Moon       [m]

    r_SC_S = r_sc_M - r_S_M;      % Sun -> spacecraft   [m]
    r_M_S  = -r_S_M;              % Sun -> Moon         [m]

    %% Frame rotations: body-fixed/component frame -> J2000
    J2000_EARTH = cspice_pxform('IAU_EARTH', 'J2000', t);
    J2000_MOON  = cspice_pxform('MOON_PA',   'J2000', t);

    %% SRP and shadowing
    % shadow_function expects:
    %   r_sun = Moon -> Sun
    %   r_sc  = Moon -> spacecraft
    % Do not pass -r_sc_M here.
    F = shadow_function(R_S, R_M, r_S_M, r_sc_M);

    % SRP() is assumed to expect Sun -> spacecraft vector.
    [aSRP, daSRP_dr, ~] = SRP(r_SC_S, eta, M, A);

    %% Earth gravity
    Cmat_E = C_mat{1};
    Smat_E = S_mat{1};

    % Spacecraft position relative to Earth, expressed in Earth-fixed frame
    r_SC_E_EF = J2000_EARTH.' * r_SC_E;

    % Earth -> Moon vector, expressed in Earth-fixed frame
    r_M_E_EF = J2000_EARTH.' * r_M_E;

    % Full Earth gravity acceleration on spacecraft
    [~, a_E_SC_EF, ddU1_EF] = potentialGradient_nm_mex( ...
        Cmat_E, Smat_E, n_maxE, r_SC_E_EF, R_E, GM_E, normalized);

    % Earth point-mass acceleration on spacecraft
    [~, a_E_SC_0_EF, ~] = potentialGradient_nm_mex( ...
        Cmat_E, Smat_E, 0, r_SC_E_EF, R_E, GM_E, normalized);

    % Earth point-mass acceleration on Moon
    [~, a_E_M_0_EF, ~] = potentialGradient_nm_mex( ...
        Cmat_E, Smat_E, 0, r_M_E_EF, R_E, GM_E, normalized);

    % Rotate Earth terms back to J2000 before combining
    a_E_SC   = J2000_EARTH * a_E_SC_EF;
    a_E_SC_0 = J2000_EARTH * a_E_SC_0_EF;
    a_E_M_0  = J2000_EARTH * a_E_M_0_EF;
    ddU1     = J2000_EARTH * ddU1_EF * J2000_EARTH.';

    % Moon-centered differential Earth acceleration
    a_Earth_2body = a_E_SC_0 - a_E_M_0;

    % Earth non-central acceleration contribution on spacecraft
    % The indirect non-central Earth term at the Moon is neglected here,
    % which is consistent with the original force-budget decomposition.
    a_Earth_oblate = a_E_SC - a_E_SC_0;

    %% Moon gravity
    Cmat_M = C_mat{2};
    Smat_M = S_mat{2};

    r_SC_M_MF = J2000_MOON.' * r_sc_M;

    [~, a_M_SC_MF, ddU2_MF] = potentialGradient_nm_mex( ...
        Cmat_M, Smat_M, n_maxM, r_SC_M_MF, R_M, GM_M, normalized);

    dU2  = J2000_MOON * a_M_SC_MF;
    ddU2 = J2000_MOON * ddU2_MF * J2000_MOON.';

    %% Sun point-mass gravity and indirect acceleration
    [a_Sun_SC, ddU3] = point_mass_accel_partials(r_SC_S, GM_S);
    a_Sun_M = -GM_S * r_M_S / norm(r_M_S)^3;

    % Moon-centered differential Sun acceleration
    a_Sun_3body = a_Sun_SC - a_Sun_M;

    %% Lunar albedo force
    % lunar_albedo_accel expects Moon -> Sun, not Sun -> Moon.
    [a_alb, da_alb_dr] = lunar_albedo_accel(r_sc_M, r_S_M, 1.3, A, M, 0.12);

    %% Relativistic corrections
    [a_rel, da_rel_dr, da_rel_dv] = relativistic_accel(r_sc_M, v_sc_M, GM_M);
    
    %% Total acceleration and partials
    a_total = dU2 + a_Earth_2body + a_Earth_oblate + a_Sun_3body + ...
              F*aSRP + F*a_alb + a_rel;

    % Approximate STM partials. The shadow-function derivative dF/dr is not
    % included because hard-shadow models are discontinuous. This is standard
    % unless using a smooth penumbra model.
    T = ddU2 + ddU1 + ddU3 + F*daSRP_dr + F*da_alb_dr + da_rel_dr;

    %% STM propagation
    PHI = reshape(x(Nx+1:Nx+Nx*Nx), [Nx, Nx]);

    J = zeros(Nx, Nx);
    J(1:3, 4:6) = eye(3);
    J(4:6, 1:3) = T;
    J(4:6, 4:6) = da_rel_dv;

    PHI_dot = J * PHI;

    %% Differential equations
    dx = [v_sc_M;
          a_total;
          PHI_dot(:)];
end

function [a, T] = point_mass_accel_partials(r, GM)
    r2 = dot(r, r);
    r1 = sqrt(r2);
    r3 = r2 * r1;
    r5 = r3 * r2;

    a = -GM * r / r3;
    T = -GM * (eye(3)/r3 - 3*(r*r.')/r5);
end

function [a_rel, da_dr, da_dv] = relativistic_accel(r, v, GM)
%RELATIVISTIC_ACCEL_PARTIALS Schwarzschild correction and partials.
    c = 299792458;   % [m/s]
    I = eye(3);

    r = r(:);
    v = v(:);

    rho = norm(r);
    v2  = dot(v, v);
    rv  = dot(r, v);

    k = GM / c^2;

    a_rel = k / rho^3 * ...
        ( (4*GM/rho - v2)*r + 4*rv*v );

    rho3 = rho^3;
    rho4 = rho^4;
    rho5 = rho^5;
    rho6 = rho^6;

    f = 4*GM/rho4 - v2/rho3;
    grad_f = (-16*GM/rho6 + 3*v2/rho5) * r;

    d_term1_dr = f*I + r*grad_f.';

    grad_rv_over_rho3 = v/rho3 - 3*rv*r/rho5;
    d_term2_dr = 4 * v * grad_rv_over_rho3.';

    da_dr = k * (d_term1_dr + d_term2_dr);

    da_dv = k * ( ...
        -2/rho3 * (r*v.') + ...
         4/rho3 * (v*r.' + rv*I) );
end

function [a_alb, da_dr] = lunar_albedo_accel(r_sc_M, r_S_M, CR, Area, mass, albedo)
%LUNAR_ALBEDO_ACCEL Simple lunar albedo acceleration and partial.
    c  = 299792458;        % [m/s]
    S0 = 1361;             % [W/m^2]
    AU = 149597870700;     % [m]
    R_M = 1738e3;          % [m]

    I = eye(3);

    r = r_sc_M(:);
    s = r_S_M(:);

    rho = norm(r);
    rhoS = norm(s);

    u_r = r / rho;
    u_s = s / rhoS;

    Sflux = S0 * (AU/rhoS)^2;
    K = CR * Area/mass * (Sflux/c) * albedo * R_M^2;

    cos_alpha = dot(u_r, u_s);
    phase_raw = 0.5 * (1 + cos_alpha);
    phase = max(0, phase_raw);

    if phase <= 0
        a_alb = zeros(3,1);
        da_dr = zeros(3,3);
        return
    end

    g = r / rho^3;
    a_alb = K * phase * g;

    dg_dr = I/rho^3 - 3*(r*r.')/rho^5;

    dcos_dr = (I - u_r*u_r.') * u_s / rho;
    dphase_dr = 0.5 * dcos_dr;

    da_dr = K * (phase * dg_dr + g * dphase_dr.');
end
