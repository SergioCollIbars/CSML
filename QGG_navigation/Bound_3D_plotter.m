clear; 
clc;
close all;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

system = "CR3BP";
[planetParams, poleParams, C_mat, S_mat, ~, ~] = ...
    load_universe(system, [0, 0], 0);
planetParams(6) = 10;

% Earth and Moon reference radius
Ref_E = planetParams(4);                                    % [m]
Ref_M = planetParams(5);                                    % [m]

% Earth & Moon distance
EarthD = planetParams(1) * planetParams(2);                 % [m]
MoonD  = (1-planetParams(1)) * planetParams(2);             % [m]

% Measurement specifications
sigmaMeas = 1E-12;                                          % [1/s^2]
R_QGG = diag([sigmaMeas,sigmaMeas,sigmaMeas,...
    sigmaMeas,sigmaMeas,sigmaMeas].^2);                     % [1/s^4]

% === 1) Define example numeric function ================================
f = @(x,y,z) computeF(x,y,z, planetParams, C_mat, S_mat, R_QGG, Ref_M, Ref_E);

% === 2) Domain/grid ====================================================
Ref = 10 * Ref_M;                        % [m]
xmin = -Ref; xmax = Ref;
ymin = -Ref; ymax = Ref;
zmin = -Ref; zmax = Ref;

N = 2000;
Nx = N; Ny = N; Nz = N;        % resolution
nLevels = 12;                  % contour levels

xv = linspace(xmin, xmax, Nx);
yv = linspace(ymin, ymax, Ny);
zv = linspace(zmin, zmax, Nz);

% === 3) Figure setup ===================================================
figure('Color','w'); hold on; box on; grid on
%axis([xmin xmax ymin ymax zmin zmax]); 
%axis vis3d tight; daspect([1 1 1])   % equal aspect, tight limits
xlabel('X [Km]'); ylabel('Y [Km]'); zlabel('Z [Km]'); view(142.8262, 10.3598)

% Helper: draw a less transparent wall surface
makeWall = @(X,Y,Z,V) surf(X,Y,Z,V, ...
    'EdgeColor','none','FaceAlpha',0.95,'Clipping','on');

% === 6) Plane z=0 (evaluated, shown on Z=zmin wall) ====================
[Xz0, Yz0] = meshgrid(xv, yv);
Vz0 = f(Xz0, Yz0, 0*Xz0);
makeWall(Xz0./1E3, Yz0./1E3, zmin*ones(size(Xz0))./1E3, Vz0);
contour3(Xz0./1E3, Yz0./1E3, zmin*ones(size(Xz0)), Vz0,  nLevels, 'k-');
disp('Finish projection onto Z = 0 plane')

% === 4) Plane x=0 (evaluated, shown on X=xmin wall) ====================
[Yx0, Zx0] = meshgrid(yv, zv);
Vx0 = f(0*Yx0, Yx0, Zx0);
makeWall(xmin*ones(size(Yx0))./1E3, Yx0./1E3, Zx0./1E3, Vx0);
contour3(xmin*ones(size(Yx0)), Yx0./1E3, Zx0./1E3, Vx0, nLevels, 'k-');
disp('Finish projection onto X = 0 plane')

% === 5) Plane y=0 (evaluated, shown on Y=ymin wall) ====================
[Xy0, Zy0] = meshgrid(xv, zv);
Vy0 = f(Xy0, 0*Xy0, Zy0);
makeWall(Xy0./1E3, ymin*ones(size(Xy0))./1E3, Zy0./1E3, Vy0);
contour3(Xy0./1E3, ymax*ones(size(Xy0))./1E3, Zy0./1E3, Vy0, nLevels, 'k-');
disp('Finish projection onto Y = 0 plane')

% ... after you compute Vx0, Vy0, Vz0 for the three walls:
Vall = [Vx0(:); Vy0(:); Vz0(:)];
Vall = Vall(~isnan(Vall));
clim([min(Vall) max(Vall)])   % << lock the color scale to the wall values
colorbar                       % (optional) show the colorbar now

% === 7) Middle point ===================================================
pc = [0,0,0];
R = Ref_M;   % radius in meters (set this to what you need)

% Sphere coordinates
[XS, YS, ZS] = sphere(50);   % resolution = 50
XS = R * XS + pc(1);         % scale to radius R, shift to center pc
YS = R * YS + pc(2);
ZS = R * ZS + pc(3);

surf(XS./1E3, YS./1E3, ZS./1E3, ...
    'FaceColor',[0.7 0.7 0.7], ...   % light gray sphere
    'EdgeColor','none');

% % load("trajectory.mat");
% % X = state_b(1, :).*planetParams(2);
% % Y = state_b(2, :).*planetParams(2);
% % Z = state_b(3, :).*planetParams(2);
% % 
% % hold on
% % hPath = plot3(X, Y, Z, '-', ...
% %     'Color', [0 0 0], ...     % solid black (won't touch colormap)
% %     'LineWidth', 2);
% % axis equal
% === 8) Title =======================================================
title('Attitude Bound in log 10 scale')


% clear kernels
cspice_kclear


%% FUNCTION
function [U] = computeF(x,y,z, planetParams, C_mat, S_mat, R_QGG, Ref_M, Ref_E)
    U = ones(length(x), length(x)) * NaN;
    for j = 1:length(x)
        for k =1:length(y)
           % position relative to the Moon
            state = [x(j, k);y(j, k);z(j, k)];
        
            % Position of the primaries (CR3BP)
            time  = 0;
            [posE, posM] = computePos_circular(time, planetParams);
        
            [ddU] = compute_gradiometer_measurements(state, posE, posM, ...
                        planetParams, C_mat, S_mat);
        
            % compute relative postion Eartn and Moon
            relPos_M = state;
            relPos_E = state - posE + posM;
        
% %             p1 = relPos_M(1); p2 = relPos_M(2); p3 = relPos_M(3);
% %         
% %             % ENU transformation
% %             lan = atan2(p2, p1);
% %             p = atan2(p3, sqrt(p1^2 + p2^2));
% %             ENU_ECEF = [-sin(lan), cos(lan), 0;...
% %                 -sin(p)*cos(lan), -sin(p)*sin(lan), cos(p);...
% %                 cos(p)*cos(lan), cos(p)*sin(lan), sin(p)];
        
            % attitude measurement partials
            T = reshape(ddU', [9, 1]); 
% %             BN = ENU_ECEF;
            BN = eye(3,3);
            [hrot] = compute_rotPartials_analy(T, BN);
            Hatt = [hrot(1, :);hrot(2, :);hrot(3, :);hrot(5, :);...
                hrot(6,:);hrot(9, :)];

% %             [Hpos] = compute_pos_partials_CR3BP(planetParams, state, ...
% %             C_mat, S_mat, posE, posM); % [1/s^2 * 1/m]
        
            if(vecnorm(relPos_M) < Ref_M || vecnorm(relPos_E) < Ref_E)
               U(j, k)   = NaN;                  % [m]
            else
                 % compute attitude Upper Bound (UB)
                f = Hatt' * (R_QGG\Hatt);
% %                 f = Hpos' * (R_QGG\Hpos);
                l = min(eig(f));
                if(rank(f) < 3)
                    BD_att = NaN;
                else
                    BD_att = ((1/l))^(1/2);
                end
               U(j, k)   = log10(BD_att);               % [rad]
            end
        end
    end
end



function [ddU] = compute_gradiometer_measurements(state, posE, posM,...
    planetParams, C_mat, S_mat)
    % extract planet paramters (non-dimensional units)
    GM_E = planetParams(8);     % [m^3 s^-2]
    GM_M = planetParams(9);     % [m^3 s^-2]
    Re_E = planetParams(4);     % [m]
    Re_M = planetParams(5);     % [m]
    
    n_max      = planetParams(6);
    normalized = planetParams(7);

    relE = state(1:3) + posM - posE;
    relM = state(1:3);
    Cmat1 = C_mat{1};
    Smat1 = S_mat{1};
    [~, ~, ddU1] = potentialGradient_nm(Cmat1, Smat1, n_max, ...
                                                relE, Re_E, GM_E, ...
                                                normalized);
    Cmat2 = C_mat{2};
    Smat2 = S_mat{2};
    [~, ~, ddU2] = potentialGradient_nm(Cmat2, Smat2, n_max, ...
                                                relM, Re_M, GM_M, ...
                                                normalized);
    
    % Gravity gradient
    ddU =  ddU1 + ddU2;
end

function [pos_earth, pos_moon] = computePos_circular(t, planetParams)
    mu = planetParams(1);
    M  = t;

    % Earth position
    pos_earth = -[mu*cos(M);mu*sin(M);0];           % [-]

    % Moon position
    pos_moon  = -[(mu-1)*cos(M);(mu-1)*sin(M); 0];  % [-]

    % bacl to meters
    pos_earth = pos_earth * planetParams(2);        % [m]
    pos_moon  = pos_moon  * planetParams(2);        % [m]
end

function [H] = compute_pos_partials_CR3BP(planetParams, x, C_mat, S_mat, posE, posM)
    eps = 1E-6;
    H = ones(6, 3) * NaN;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = x + Ar./2;   % [ACI]
        rneg = x - Ar./2;   % [ACI]

        [ddU_pos] = compute_gradiometer_measurements(rpos, posE, posM, ...
                        planetParams, C_mat, S_mat);

        [ddU_neg] = compute_gradiometer_measurements(rneg, posE, posM, ...
                        planetParams, C_mat, S_mat);

        Ht = (ddU_pos - ddU_neg)./(vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end
