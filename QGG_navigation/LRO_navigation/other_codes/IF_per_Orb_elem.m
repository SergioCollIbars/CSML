clear;
clc;
close all;
addpath('simplified_functions/');
addpath('data/');
format long g;
cspice_furnsh(...
    '/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')

%% INFORMATION MATRIX PER INCLINATION
% Compute the information matrix per inclination for a Lunar orbit
% Date: 10/10/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Body constants (Moon defaults; change as needed) ----
R   = 1737.4E3;                 % [m] 
mu  = 4.902800118E3;            % [km^3/s^2] 
[planetParams, Cmat_true, Smat_true] = load_universe();

% time
utc_start = '2025-05-20 00:00:00';
utc_stop  = '2025-05-20 02:00:00';
N         = 2000;           % number of samples
et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

% filter
fCut = 1E-4;
Fs = 1 / (et(2) - et(1));
d = designfilt('highpassiir', ...
    'FilterOrder', 12, ...
    'HalfPowerFrequency', fCut, ...
    'SampleRate', Fs);

% ---- Choose ONE set of orbital elements (radians) ----
a     = R/1E3 + 90;              % [km]
e     = 0;                       % set 0 for circular test
i_range = 0:1:90;                % inclination [deg]
Omega = deg2rad(0);              % RAAN [rad]
omega = deg2rad(0);              % arg. perigee [rad]
M     = deg2rad(0);              % mean anomaly [rad]

IF_inc_pos = ones(length(i_range), length(et)) * NaN;
IF_inc_vel = ones(length(i_range), length(et)) * NaN;

%% compute inclination = 0
tic;
    incIdx = 1;
    disp('Computing inclination ...  ' +string(i_range(incIdx)));
    i = deg2rad(i_range(incIdx));
    
    % Compute initial conditions. X0
    [rECI, vECI] = coe2rv(a,e,i,Omega,omega,M,mu);
    X0 = [rECI;vECI].*1E3;  % [m]
    
    % integrate trajectory
    options = odeset('RelTol',1E-13,'AbsTol',1E-13);
    STM0 = reshape(eye(6,6), [36, 1]);
    [t, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
        Cmat_true, Smat_true), et, [X0; STM0], options);
    
    % compute RTN matrix
    BN_mat = ones(3*length(t), 3)*NaN;
    for k = 1:length(t)
        maxInd = 3 * k; minInd = maxInd - 2;
        [NB] = RTN2ECI(state(k, 1:3)', state(k, 4:6)');
        BN_mat(minInd:maxInd, :) = NB';
    
    end
    
    [sigma_time] = compute_IF(t, state, Cmat_true, ...
    Smat_true, BN_mat, d);
    
    IF_inc_pos(incIdx, :) = sum(sigma_time(1:3, :)); 
    IF_inc_vel(incIdx, :) = sum(sigma_time(4:6, :));
dt = toc;
expected_time = dt * 90 / 60;
disp('      Expected computation time ... ' + string(expected_time) + ' min');

%% compute rest of coefficients
 for incIdx = 2:length(i_range)
     disp('Computing inclination ...  ' +string(i_range(incIdx)));
     i = deg2rad(i_range(incIdx));

    % Compute initial conditions. X0
    [rECI, vECI] = coe2rv(a,e,i,Omega,omega,M,mu);
    X0 = [rECI;vECI].*1E3;  % [m]
    
    % integrate trajectory
    options = odeset('RelTol',1E-13,'AbsTol',1E-13);
    STM0 = reshape(eye(6,6), [36, 1]);
    [t, state] = ode113(@(t, x) EOM_LRO_EPHEM(t, x, planetParams, ...
        Cmat_true, Smat_true), et, [X0; STM0], options);
    
% %     figure()
% %     plot3(state(:, 1), state(:, 2), state(:, 3), 'LineWidth', 2);
% %     axis equal;
% %     hold on;
% % 
% %     % Create a sphere
% %     [Xs, Ys, Zs] = sphere(100);      % resolution 100x100
% %     surf(R*Xs, R*Ys, R*Zs, ...
% %         'FaceAlpha', 1, ...          % transparency
% %         'EdgeColor', 'none', ...
% %         'FaceColor', [0.8 0.8 0.8]); % light gray
% %     
% %     axis equal;
% %     xlabel('X [km]')
% %     ylabel('Y [km]')
% %     zlabel('Z [km]')
% %     title('LRO Trajectory around the Moon. J2000 frame')
    
    % compute RTN matrix
    BN_mat = ones(3*length(t), 3)*NaN;
    for k = 1:length(t)
        maxInd = 3 * k; minInd = maxInd - 2;
        [NB] = RTN2ECI(state(k, 1:3)', state(k, 4:6)');
        BN_mat(minInd:maxInd, :) = NB';
    end

    [sigma_time] = compute_IF(t, state, Cmat_true, ...
    Smat_true, BN_mat, d);

    IF_inc_pos(incIdx, :) = sum(sigma_time(1:3, :)); 
    IF_inc_vel(incIdx, :) = sum(sigma_time(4:6, :));
 end
 
figure()
[x, y] = meshgrid(t, i_range);
contourf(x, y, real(IF_inc_pos))
xlabel('time')
ylabel('inclination [deg]')
colorbar();
title('Position trace uncertainty')

figure()
contourf(x, y, real(IF_inc_vel))
xlabel('time')
ylabel('inclination [deg]')
colorbar();
title('Velocity trace uncertainty')

% uncertainty trace at final time
figure()
subplot(1, 2, 1)
plot(i_range, 3.*sum(IF_inc_pos(:, end)'), 'LineWidth', 2)
xlabel('inclination [deg]');
ylabel('[m]');
title('Position trace 3\sigma')
grid on;

subplot(1, 2, 2)
plot(i_range, 3.*sum(IF_inc_vel(:, end)'), 'LineWidth', 2)
xlabel('inclination [deg]');
ylabel('[m/s]');
title('Velocity trace 3\sigma')
grid on;


%% FUNCTIONS
function [rI, vI] = coe2rv(a,e,i,Omega,omega,M,mu)
    % Elements -> ECI r,v (scalar)
    E = keplerE(M,e);
    f = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
    p = a*(1 - e^2);
    r = p/(1 + e*cos(f));
    r_pf = [r*cos(f); r*sin(f); 0];
    v_pf = sqrt(mu/p)*[-sin(f); e+cos(f); 0];
    
    Q = rotz(Omega)*rotx(i)*rotz(omega);
    rI = Q*r_pf; vI = Q*v_pf;
end

function E = keplerE(M,e)
    % Solve Kepler's equation M = E - e sin E (scalar)
    M = mod(M, 2*pi);
    if e < 0.8, E = M; else, E = pi; end
    for k=1:50
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f/fp;
        E  = E + dE;
        if abs(dE) < 1e-14, break; end
    end
end

function R = rotx(a)
    ca=cos(a); sa=sin(a);
    R = [1 0 0; 0 ca -sa; 0 sa ca];
    end

function R = rotz(a)
    ca=cos(a); sa=sin(a);
    R = [ca -sa 0; sa ca 0; 0 0 1];
end

function [H] = compute_partials_gradiometer(t, k, state, Cmat_true, Smat_true, BN_mat)
    % rotation matrix
    frame_from = 'MOON_PA';
    frame_to   = 'J2000';
    J2000_MOON = cspice_pxform(frame_from, frame_to, t(k));

    % Moon parameters
    [GM2] = cspice_bodvrd('MOON', 'GM', 1);    % Get GM for the Sun [km^3/s^2]
    [Re2] = cspice_bodvrd('MOON', 'RADII', 3); % get Moon Radius [Km]
    GM2 = GM2*1E9;                             % [m^3 s^-2]
    Re2 = Re2*1E3;                             % [m]

    % state in the S/C frame
    maxInd = 3 * k; minInd = maxInd - 2;
    BN = BN_mat(minInd:maxInd, :);
    x = BN * state(k, 1:3)';

    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [S/C FRAME]
        rneg = x - Ar./2;   % [S/C FRAME]

        % rotation from inertial S/C frame to Moon frame
        Rot = BN * J2000_MOON;


        Cmat2 = Cmat_true{2};
        Smat2 = Smat_true{2};
        [~, ~, ddU_pos_ECEF] = potentialGradient_nm(Cmat2, Smat2, 60, ...
                                                Rot'*rpos, Re2(1), GM2, ...
                                                1);
        [~, ~, ddU_neg_ECEF] = potentialGradient_nm(Cmat2, Smat2, 60, ...
                                                Rot'*rneg, Re2(1), GM2, ...
                                                1);

        [~, ~, ddU0_pos_ECEF] = potentialGradient_nm(Cmat2, Smat2, 0, ...
                                                Rot'*rpos, Re2(1), GM2, ...
                                                1);
        [~, ~, ddU0_neg_ECEF] = potentialGradient_nm(Cmat2, Smat2, 0, ...
                                                Rot'*rneg, Re2(1), GM2, ...
                                                1);
        
        ddU_pos = Rot * (ddU_pos_ECEF - ddU0_pos_ECEF) * Rot';
        ddU_neg = Rot * (ddU_neg_ECEF - ddU0_neg_ECEF) * Rot';

% %         ddU_pos = Rot * (ddU_pos_ECEF) * Rot';
% %         ddU_neg = Rot * (ddU_neg_ECEF) * Rot';
        Ht      = (ddU_pos - ddU_neg)./(vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end

function [sigma_time] = compute_IF(t, state, Cmat_true, ...
    Smat_true, BN_mat, d)

    % measurement weight
    Rinv = 1/(1E-24) * eye(6); 

    % STM at different times
    STM = state(:, 7:end);

    IF0 = 0;
    H_stack = zeros(length(t), 6, 6);
    % compute information matrix
    for k = 1:length(t)
        % STM from current to t0
        PHI = reshape(STM(k, :), [6, 6]);
    
       % compute measurements at current time
       [H_pos] = compute_partials_gradiometer(t, k, state, ...
           Cmat_true, Smat_true, BN_mat);
       H_k_0 = [H_pos, zeros(6, 3)] * PHI;
       
       H_stack(k, :, :) = H_k_0;

       % information increment
       delta_IF = H_k_0' * Rinv * H_k_0;
       
       % update information
       IF0 = IF0 + delta_IF;
    end

    % select measurement & filter in frequency
    Nm = 6;
    H = zeros(Nm * length(t), 6);
    for k = 1:Nm
        h = squeeze(H_stack(:, k, :));
        Hf = filtfilt(d, h);   % Hf is [N x n]
        maxInd = k * length(t); minInd = maxInd - (length(t)-1);
        H(minInd:maxInd, :) = Hf;
    end

    C = repmat({Rinv}, 1, length(t));   % Create 1xN cell array of A
    M = blkdiag(C{:});       % Expand the cell contents into blkdiag
    IF0 = H' * M * H;
   
    % compute uncertainty over time
    P0 = inv(IF0);
    sigma_time = ones(6, length(t)) * NaN;
    for k = 1:length(t)
        PHI = reshape(STM(k, :), [6, 6]);
        Pk = PHI * P0 * PHI';
        sigmak = sqrt(diag(Pk));

        % save uncertainty at current time
        sigma_time(:, k) = sigmak';
    end
% %     PHI = reshape(STM(end, :), [6, 6]);
% %     IF = inv(PHI)' * IF0 * inv(PHI);
% %     sigma = sqrt(inv(IF));
% % 
% %     % compute stats
% %     IF_pos = [sigma(1,1);sigma(2,2);sigma(3,3)];
% %     IF_vel = [sigma(4,4);sigma(5,5);sigma(6,6)];
end