clear;
clc;
close all;
addpath("functions/")
set(0,'defaultAxesFontSize',16);
%%         GRADIOMETER POSITION PARTIALS SENSITIVITY TEST
% Description: Test functions to compute gradiometer partials to position
%   compute partials in the inertial frame.
%   Also test gradiometer invariants partials to position in the iniertial
%   frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% planet mass parameters
GM_E = 3.986004E14;     % [m^3/s^2]
GM_M = 4.9048695E14;    % [m^3/s^2]
GM_B = 5.2;             % [m^3/s^2]

% radius parameters
R_E =  6.3781E6;        % [m]
R_M =  1740E3;          % [m]
R_B = 250;              % [m]

% SH perturbations (not important. nMax = 0)
Cmat = [1, 0, 0, 0;...
     0, 0, 0, 0;...
    -1.08E-3, 0, 1.57E-6, 0;...
     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 

Smat = [0, 0, 0, 0;...
     0, 0, 0, 0;...
     0, 0, -9.03E-7, 0;...
     0, 2.68E-7, -2.12E-7, 1.98E-7]; 

% select planet parameters
GM = GM_M;
R = R_M;

% position vector
x = 3000;
y = 500;
z = 750;

% position deviation
N = 1000;
delta = linspace(1, 1000, N);
ths = 1E-16;
xlb = ["\Gamma_{xx}", "\Gamma_{xy}", "\Gamma_{xz}", "\Gamma_{yy}", ...
    "\Gamma_{yz}", "\Gamma_{zz}"];
% output matrices
error = ones(6, N) * NaN;
for k = 1:3
    w = zeros(3, 1);
    for j = 1:N
        w(k) = 1;
        d = delta(j) * w;

        % compute grad meas at nominal
        r = [x;y;z];
        [~, ddU_nom] = gradiometer_meas(0, r, Cmat, Smat, ...
        [0,0,0,0], [GM, R, 0, 0]);
        ddU_nom = reshape(ddU_nom, [3,3]);
    
        % compute grad meas at deviation
        r = [x;y;z] + d;
        [~, ddU_dev] = gradiometer_meas(0, r, Cmat, Smat, ...
        [0,0,0,0], [GM, R, 0, 0]);
        ddU_dev = reshape(ddU_dev, [3,3]);
    
        % compute partials
        [ddU_dxyz] = compute_grad_posPartials(GM, x, y, z);
    
        % compute error
        if(abs(ddU_dev(1,1) - ddU_nom(1,1)) < ths)
            error(1, j) = 0;
        else
            error(1, j) = (ddU_dev(1,1) - ddU_nom(1,1)) - ...
                ddU_dxyz(1, :) * d;
            error(1, j) =  error(1, j) / (ddU_dev(1,1) - ddU_nom(1,1));
        end
        
        if(abs(ddU_dev(1,2) - ddU_nom(1,2)) < ths)
            error(2, j) = 0;
        else
            error(2, j) = (ddU_dev(1,2) - ddU_nom(1,2)) - ...
                ddU_dxyz(2, :) * d;
            error(2, j) =  error(2, j) / (ddU_dev(1,2) - ddU_nom(1,2));
        end
    
        if(abs(ddU_dev(1,3) - ddU_nom(1,3)) < ths)
            error(3, j) = 0;
        else
            error(3, j) = (ddU_dev(1,3) - ddU_nom(1,3)) - ...
                ddU_dxyz(3, :) * d;
            error(3, j) =  error(3, j) / (ddU_dev(1,3) - ddU_nom(1,3));
        end
    
        if(abs(ddU_dev(2,2) - ddU_nom(2,2)) < ths)
            error(4, j) = 0;
        else
            error(4, j) = (ddU_dev(2,2) - ddU_nom(2,2)) - ...
                ddU_dxyz(4, :) * d;
            error(4, j) =  error(4, j) / (ddU_dev(2,2) - ddU_nom(2,2));
        end
    
        if(abs(ddU_dev(2,3) - ddU_nom(2,3)) < ths)
            error(5, j) = 0;
        else
            error(5, j) = (ddU_dev(2,3) - ddU_nom(2,3)) - ...
                ddU_dxyz(5, :) * d;
            error(5, j) =  error(5, j) / (ddU_dev(2,3) - ddU_nom(2,3));
        end
    
        if(abs(ddU_dev(3,3) - ddU_nom(3,3)) < ths)
            error(6, j) = 0;
        else
            error(6, j) = (ddU_dev(3,3) - ddU_nom(3,3)) - ...
                ddU_dxyz(6, :) * d;
            error(6, j) =  error(6, j) / (ddU_dev(3,3) - ddU_nom(3,3));
        end
    end
    figure()
    for h = 1:6
        subplot(2, 3, h)
        plot(delta, error(h, :), 'LineWidth', 2, 'Color', 'r')
        xlabel('\delta x_i [m]');
        ylabel('% error ' + xlb(h))
    end
    sgtitle('1st order approx. Error for \delta x_' + string(k));
end
