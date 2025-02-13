clc;
clear;
close all;
format long g;

%%          ATTITUDE ERROR ANALYSIS
%   Description: Analysis on how the attitude influence the measurements.

% IMPORTS
addpath('/Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/');
addpath('/Users/sergiocollibars/Desktop/CSML/codes/QGG/src/');
addpath('functions/');

% read position in B coords
filename = 'orbitData.txt';
posData = importdata(filename).data;
r_B = [posData(:, 5), posData(:, 6), posData(:, 7)]';
v_B = [posData(:, 11), posData(:, 12), posData(:, 13)]';

% read position in ACI coords
r_ACI = [posData(:, 2), posData(:, 3), posData(:, 4)]';

% read measurements in B coords
filename = 'accData.txt';
accData = importdata(filename).data;
accT = [accData(:, 2), accData(:, 5), accData(:, 8), ...
    accData(:, 3), accData(:, 6), accData(:, 9), ...
    accData(:, 4), accData(:, 7), accData(:, 10)]';

% read angular velocity and angular acceleration
filename = 'attitudeData.txt';
attData = importdata(filename).data;
omega = [attData(:, 11), attData(:, 12), attData(:, 13)]';
Omega = [attData(:, 14), attData(:, 15), attData(:, 16)]';

% read time 
TIME = accData(:, 1);

% Rotation matrices
filename = 'ECI2Body.txt';
BN = importdata(filename).data;

filename = 'N2ACAF.txt';
ACAF_N = importdata(filename).data;

% Harmonics values
Cnm = [1 0 0 0 0 0 0;...
                    0 0 0 0 0 0 0;...
                    -0.0391557863539988 -2.96928723209235e-06 0.00375640654748657 0 0 0 0;...
                    0.0148427177700986 0.00167095097673949 3.80845003468165e-05 0.000371755938456641 0 0 0;...
                    0.0307494000000000 0.000413625917950024 -0.000490123739988179 -6.43092753252371e-05 4.51227856599555e-05 0 0;...
                    -0.00456599734888228 0.000441961635589938 8.84264903209225e-05 1.40396087936725e-05 1.07458402471302e-05 1.19886311472876e-06 0;...
                    -0.00896736657720649 0.000905916675449317 9.05780703119059e-05 2.77025499573633e-05 -1.92680576794137e-06 9.97533032070693e-08 1.67034838314692e-07];

Snm = [0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            0 0 -2.54325906400954e-05 0 0 0 0;...
            0 0.000992000134408593 6.53000000000000e-05 -0.00100797120329237 0 0 0;
            0 0.000634013000392474 0.000108054642426876 0.000102400000000000 0.00291093983173820 0 0;...
            0 -1.75754943031483e-05 -7.41878397813858e-05 4.79413750994751e-06 -0.000503800000000000 0.000448812426298560 0;...
            0 -1.35941163743731e-06 -9.69889209840526e-06 -7.55736000569855e-06 -1.59676058718751e-06 -0.000599300000000000 -3.93397896234722e-05];


% Bennu parameters
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
T = 2 * pi / n;

% Position error. B
deltaR_B = [10; 0; 0];     % [m]

% simulation index
Nt_max = 77760;

% matrix definition
rP_B = zeros(3, Nt_max);
rP_ACI = zeros(3, Nt_max);
rP_ACAF = zeros(3, Nt_max);

Err = zeros(9, Nt_max);
Err_ENU = zeros(9, Nt_max);

% loop
for Nt = 1:Nt_max
    % Perturbed position error. ACI & ACAF
    up = Nt * 3;
    down = up - 2;
    
    rP_B(:, Nt) = r_B(:, Nt) + deltaR_B;
    rP_ACI(:, Nt) = BN(down:up, :)' * rP_B(:, Nt);
    rP_ACAF(:, Nt) = ACAF_N(down:up, :) * rP_ACI(:, Nt);
    
    % compute angular velocity. B
    [omega2M, OmegaM] = angularMatrix(omega(:, Nt), Omega(:, Nt));
    
    % compute gravity tensor
    [~, dUp_ACAF, ddUp_ACAF] = potentialGradient_nm(Cnm, Snm, 6, ...
                                                    rP_ACAF(:, Nt), Re, GM, 0);
    
    ddUp_ACI = ACAF_N(down:up, :)' * ddUp_ACAF * ACAF_N(down:up, :);
    ddUp_B = BN(down:up, :) * ddUp_ACI * BN(down:up, :)';
    
    % compute measurement
    accTp = -ddUp_B + OmegaM + omega2M; 
    accTp = reshape(accTp, [9, 1]);
    
    % compute error. B components (RTN)
    Err(:, Nt) = (accT(:, Nt) - accTp);

    % longitude latitude.
    lambda = atan2(rP_ACAF(2, Nt), rP_ACAF(1, Nt)); 
    phi = atan2(rP_ACAF(3, Nt), sqrt(rP_ACAF(2, Nt)^2 + rP_ACAF(1, Nt)^2));

    ENU_ACAF = [-sin(lambda), cos(lambda), 0;...
            -cos(lambda)*sin(phi), -sin(lambda)*sin(phi), cos(phi);...
            cos(lambda)*cos(phi), sin(lambda)*cos(phi), sin(phi)];
    B_ENU = BN(down:up, :) * ACAF_N(down:up, :)' * ENU_ACAF';
    
    err_ENU = B_ENU' * reshape(Err(:, Nt), [3,3]) * B_ENU;
    Err_ENU(:, Nt) = reshape(err_ENU, [9, 1]);

end

% orbit plot
plotOrbit(r_ACI, rP_ACI, Nt_max);

% plot error over time
ttitle = 'Error along orbit. RTN components';
ytitle = ["R", "RT", "RN", "T", "TN", "N"];
plotT(Err, Nt_max, TIME, ttitle, ytitle);

ttitle = 'Error along orbit. ENU components';
ytitle = ["E", "EN", "EU", "N", "NU", "U"];
plotT(Err_ENU, Nt_max, TIME, ttitle, ytitle);

% plot tensor
ttitle = 'Sensor truth measurements';
ytitle = ["R", "RT", "RN", "T", "TN", "N"];
plotT(accT(:,1:Nt_max), Nt_max, TIME, ttitle, ytitle);


%%                      FUNCTIONS
function [omega, Omega] = computeAng(ro, vo)
    theta_dot = cross(ro, vo) / (vecnorm(ro)^2);
    theta_dot_u = theta_dot / vecnorm(theta_dot); 

    h = cross(ro, vo);
    Omega  = -2 * vecnorm(h) / (vecnorm(ro)^4) * ...
        dot(ro, vo) * theta_dot_u;
    omega = theta_dot;
end

function plotOrbit(r_ACI, rP_ACI, Nt)
    % plot 3D orbit
    figure();
    plot3(r_ACI(1,1:Nt), r_ACI(2,1:Nt), r_ACI(3,1:Nt), 'MarkerSize', 10);
    hold on;
    plot3(rP_ACI(1, 1:Nt), rP_ACI(2, 1:Nt), rP_ACI(3, 1:Nt), 'MarkerSize', 10);
    legend('Nominal', 'Perturbed')
    title('Orbit in the inertial frame {I, J, K}');
    grid on;
    
    % plot axis
    mAxis = max(max(r_ACI));
    %axis([0 max(ri(1,:)) 0 max(ri(2,:)) 0 max(ri(3,:))])
    axis([0 mAxis 0 mAxis 0 mAxis])
    hold all;
    quiver3(0,0,-max(0),0,0,max(zlim),'r','LineWidth',1)
    quiver3(0,-max(0),0,0,max(ylim),0,'r','LineWidth',1)
    quiver3(-max(0),0,0,max(xlim),0,0,'r','LineWidth',1)
    text(0,0,max(zlim),'K','Color','r')
    text(0,max(ylim),0,'J','Color','r')
    text(max(xlim),0,0,'I','Color','r')
    
    scale =  450;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 
    
    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end

function plotT(Err, Nt_max, TIME, ttitle, ytitle)
    % time vector [h]
    t = TIME(1:Nt_max)./ 87132.103;

    figure();

    subplot(2, 3, 1);
    plot(t, Err(1, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(1))
    grid on;

    subplot(2, 3, 2);
    plot(t, Err(2, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(2))
    grid on;

    subplot(2, 3, 3);
    plot(t, Err(3, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(3))
    grid on;

    subplot(2, 3, 4);
    plot(t, Err(5, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(4))
    grid on;

    subplot(2, 3, 5);
    plot(t, Err(6, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(5))
    grid on;

    subplot(2, 3, 6);
    plot(t, Err(9, :))
    xlabel("Orbital rev, T")
    ylabel(ytitle(6))
    grid on;

    sgtitle(ttitle)
end