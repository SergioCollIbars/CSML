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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       CRAMER-RAO UPPER BOUND                            %
% Description: Evaluate the CRUB in the cislunar space using the CR3BP    %
% model.                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% create X Y domain
Ref   = planetParams(2);                                    % [m]
% % Ref = Ref_M;
N     = 50;
valsX = [-1.5*Ref, 1.5.*Ref, N];             % [m]
valsY = [-1.5*Ref, 1.5*Ref, N];                             % [m]

IF_QGG_pos   = ones(valsX(3), valsY(3)) * NaN;
BD_QGG_pos   = ones(valsX(3), valsY(3)) * NaN;
IF_QGG_att   = ones(valsX(3), valsY(3)) * NaN;
BD_QGG_att   = ones(valsX(3), valsY(3)) * NaN;

Nx = valsX(3); Ny = valsY(3);

[X, Y] = meshgrid(linspace(valsX(1), valsX(2), valsX(3)), ...
                    linspace(valsY(1), valsY(2), valsY(3)));

% Measurement specifications
sigmaMeas = 1E-12;                                          % [1/s^2]
R_QGG = diag([sigmaMeas,sigmaMeas,sigmaMeas,...
    sigmaMeas,sigmaMeas,sigmaMeas].^2);                     % [1/s^4]

% Compute Bounds over domain 
for j = 1:Nx    % x domain
    disp('Computing progress ...' + string(j/Nx * 100) + ' %');
    for i =1:Ny % y domain
        % absolute position
        state = [X(i,j), Y(i, j), 0];   % [m]

        % Position of the primaries (CR3BP)
        time  = 0;
        [posE, posM] = computePos_circular(time, planetParams);

       [ddU] = compute_gradiometer_measurements(state, posE, posM, ...
                    planetParams, C_mat, S_mat);

        % compute relative postion Eartn and Moon
        relPos_M = state' - posM;
        relPos_E = state' - posE;

        x = relPos_M(1); y = relPos_M(2); z = relPos_M(3);
        
        % position measurement partials
        [Hpos] = compute_pos_partials_CR3BP(planetParams, state', ...
            C_mat, S_mat, posE, posM); % [1/s^2 * 1/m]

        % ENU transformation
        lan = atan2(y, x);
        p = atan2(z, sqrt(x^2 + y^2));
        ENU_ECEF = [-sin(lan), cos(lan), 0;...
            -sin(p)*cos(lan), -sin(p)*sin(lan), cos(p);...
            cos(p)*cos(lan), cos(p)*sin(lan), sin(p)];

        % attitude measurement partials
        T = reshape(ddU', [9, 1]); 
        BN = ENU_ECEF;
% %         BN = eye(3,3);
        [hrot] = compute_rotPartials_analy(T, BN);
        Hatt = [hrot(1, :);hrot(2, :);hrot(3, :);hrot(5, :);...
            hrot(6,:);hrot(9, :)];
        
        % compute position Upper Bound (UB)
        f = Hpos' * (R_QGG\Hpos);
        l = min(eig(f));
        if(rank(f) < 3)
            BD_pos  = NaN;
            IF_pos  = NaN;
        else
            BD_pos  = ((1/l))^(1/2);
            IF_pos  = det(f);
        end

        % compute attitude Upper Bound (UB)
        f = Hatt' * (R_QGG\Hatt);
        l = min(eig(f));
        if(rank(f) < 3)
            BD_att = NaN;
            IF_att  = NaN;
        else
            BD_att = ((1/l))^(1/2);
            IF_att  = det(f);
        end

        if(vecnorm(relPos_M) < Ref_M || vecnorm(relPos_E) < Ref_E)
            BD_QGG_pos(i, j)   = NaN;                  % [m]
            IF_QGG_pos(i, j)   = NaN;                  % [m^-3]
            BD_QGG_att(i, j)   = NaN;                  % [m]
            IF_QGG_att(i, j)   = NaN;                  % [rad^-3]
        else
            BD_QGG_pos(i, j)   = BD_pos;               % [m]
            IF_QGG_pos(i, j)   = IF_pos;               % [m^-3]
            BD_QGG_att(i, j)   = BD_att;               % [rad]
            IF_QGG_att(i, j)   = IF_att;               % [rad^-3]
        end
    end
end


% Plot Bounds (Z = 0 plane)
Y_slice = squeeze(Y(:, :));
X_slice = squeeze(X(:, :));
IF_QGG_slice = squeeze(IF_QGG_pos(:, :));
BD_QGG_slice = squeeze(BD_QGG_pos(:, :));

tt = 'Position Bound [m]';
plot_contour(X_slice./1E3, Y_slice./1E3, BD_QGG_slice, tt)

tt = 'Position Information [1 / m^3]';
plot_contour(X_slice./1E3, Y_slice./1E3, IF_QGG_slice, tt)

IF_QGG_slice = squeeze(IF_QGG_att(:, :));
BD_QGG_slice = squeeze(BD_QGG_att(:, :));

tt = 'Attitude Bound [rad]';
plot_contour((X_slice - MoonD)./Ref_M, Y_slice./Ref_M, real(BD_QGG_slice), tt)

tt = 'Attitude Information [1 / rad^3]';
plot_contour(X_slice./1E3, Y_slice./1E3, IF_QGG_slice, tt)

% clear kernels
cspice_kclear

%% FUNCTIONS
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

        [ddU_pos] = compute_gradiometer_measurements(rpos', posE, posM, ...
                        planetParams, C_mat, S_mat);

        [ddU_neg] = compute_gradiometer_measurements(rneg', posE, posM, ...
                        planetParams, C_mat, S_mat);

        Ht = (ddU_pos - ddU_neg)./(vecnorm(rpos-rneg));
        
        H(:, j) = [Ht(1,1);Ht(1,2);Ht(1,3);Ht(2,2);Ht(2,3);Ht(3,3)];
    end
end


function [Hpos] = compute_att_partials_CR3BP(planetParams, x, C_mat, S_mat, posE, posM) % WARNING: not using this function. Only for testing!
    eps = 1E-6;
    Hpos = ones(9, 3) * NaN;
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;

        Atpos = At./2;
        Atneg = - At./2; 

        [Rpos] = rotationMatrix(Atpos(1), Atpos(2), Atpos(3), [3, 2, 1]);
        [Rneg] = rotationMatrix(Atneg(1), Atneg(2), Atneg(3), [3, 2, 1]);

        [ddU] = compute_gradiometer_measurements(x, posE, posM, ...
                        planetParams, C_mat, S_mat); % [ACI frame]

        ddU_pos = Rpos * ddU * Rpos';
        ddU_neg = Rneg * ddU * Rneg';


        H = (ddU_pos - ddU_neg)./(vecnorm(Atpos-Atneg));
        
        Hpos(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
            H(3,1); H(3,2);H(3,3)];
    end
end

function [ddU] = compute_gradiometer_measurements(state, posE, posM, ...
    planetParams, C_mat, S_mat)
    % extract planet paramters (non-dimensional units)
    GM_E = planetParams(8);     % [m^3 s^-2]
    GM_M = planetParams(9);     % [m^3 s^-2]
    Re_E = planetParams(4);     % [m]
    Re_M = planetParams(5);     % [m]
    
    n_max      = planetParams(6);
    normalized = planetParams(7);

    relE = state(1:3)' - posE;
    relM = state(1:3)' - posM;
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

function [] = plot_contour(X, Y, value, tt)
    mu    = 0.0121505843958292; 
    scale = 384399;  % [km]

    figure;
    contourf(X, Y, log10(value), 20, 'LineStyle', 'none');
% %     xlabel('X [km]');
% %     ylabel('Y [km]');
     xlabel('X [R Moon]');
    ylabel('Y [R Moon]');
    colormap('parula')
    title(tt);
    axis equal
    colorbar
    hold all;
% %     plot(-mu*scale,0, "o",'MarkerFaceColor','r', ...
% %         'MarkerEdgeColor', 'r')
% %     plot((1-mu)*scale,0,"pentagram",'MarkerFaceColor','r',...
% %         'MarkerEdgeColor', 'r')
% %     legend('','Earth', 'Moon')
end