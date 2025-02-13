close all;
clear;
clc;

%%              HARMONIC INFLUENCE IN LON LAT

% IMPORTS
addpath('functions/')
addpath('dataFiles/')

% load data
LON = load('lon.mat').lon;
LAT = load("lat.mat").lat;
N = 300;

% longitude meshgrid
lat = linspace(-pi/2, pi/2, N);
lon = linspace(-pi, pi, N);
[Px, Py] = meshgrid(lon, lat);

% Planet parameters
Re = 246;                           % Planet radius [m]   6355E3 (polar radius)
GM = 5.2;                           % gravity param [m^3 / s^2]
H = 1000;                           % SC altitude [m]

n_max = 5;                          % spherical harmonic order
 Nc = -2;
for j = 1:n_max+1
    Nc = Nc + j;
end

% Hamonics coefficients
C_mat = [1 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            -0.017511 -0.0000023 0.0058194 0 0 0 0;...
            0.00561002 0.0015471 0.0001115 0.0026660 0 0 0;...
            0.0102498 0.0004360 -0.0021919 -0.0010761 0.0021356 0 0;...
            -0.0013767 0.0005161 0.0005464 0.0004250 0.0013801 0.0004869 0;...
            -0.0024871 0.0011514 0.0007281 0.0013361 -0.0005090 0.0001236 0.00071695];

S_mat = [0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0;...
            0 0 -0.0000197 0 0 0 0;...
            0 0.0015368 0.0000653 -0.0009332 0 0 0;
            0 0.0018562 0.0007749 0.0001024 0.0030684 0 0;...
            0 -0.0000786 -0.0012414 0.0002269 -0.0005038 0.0005241 0;...
            0 -0.0000084 -0.0002936 -0.0009706 -0.0006485 -0.0005993 -0.0000500];
% Harmonic 
U = zeros(N, N, n_max+1);
R = H + Re;
disp('n = ' + string(n_max));
index = Nc - n_max;
for nn = 1:n_max+1
    disp('m = ' + string(nn-1))
    for j = 1:N
        for i =1:N
            % get current position. ACAF frame
            x = R * cos(Py(i, j)) * cos(Px(i, j));
            y = R * cos(Py(i, j)) * sin(Px(i, j));
            z = R * sin(Py(i, j));

            r_ACAF = [x;y;z];
            r_ACAF_u = r_ACAF./vecnorm(r_ACAF);

            % compute visibility matrix. ACAF frame
            [H] = potentialGradient_Cnm(C_mat, S_mat, n_max, ...
                                                r_ACAF, Re, GM, eye(3));
            A = reshape(H(:, index), [3, 3]);
            
            % gravity potential in the radius direction
            U(i, j, nn) = vecnorm(A * r_ACAF_u);
        end
    end
    index = index + 1;
end

% plot resutls
figure();
[p,n]=numSubplots(n_max+1);

for j = 1:n_max+1
    subplot(p(1), p(2), j);
    h = pcolor(rad2deg(Px), rad2deg(Py), U(:, :, j));
    set(h, 'EdgeColor', 'none');
    view(2);
    colorbar();
    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    
    hold on;
    plot(rad2deg(LON), rad2deg(LAT), 'r*', 'MarkerSize',0.1);
    title('U n = ' + string(n_max) + ' , ' + 'm = ' + string(j -1));
end

%% FUNCTIONS
function [p,n]=numSubplots(n)
    % function [p,n]=numSubplots(n)
    %
    % Purpose
    % Calculate how many rows and columns of sub-plots are needed to
    % neatly display n subplots. 
    %
    % Inputs
    % n - the desired number of subplots.     
    %  
    % Outputs
    % p - a vector length 2 defining the number of rows and number of
    %     columns required to show n plots.     
    % [ n - the current number of subplots. This output is used only by
    %       this function for a recursive call.]
    %
    %
    %
    % Example: neatly lay out 13 sub-plots
    % >> p=numSubplots(13)
    % p = 
    %     3   5
    % for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end 
    %
    %
    % Rob Campbell - January 2010
       
        
    while isprime(n) & n>4, 
        n=n+1;
    end
    p=factor(n);
    if length(p)==1
        p=[1,p];
        return
    end
    while length(p)>2
        if length(p)>=4
            p(1)=p(1)*p(end-1);
            p(2)=p(2)*p(end);
            p(end-1:end)=[];
        else
            p(1)=p(1)*p(2);
            p(2)=[];
        end    
        p=sort(p);
    end
    %Reformat if the column/row ratio is too large: we want a roughly
    %square design 
    while p(2)/p(1)>2.5
        N=n+1;
        [p,n]=numSubplots(N); %Recursive!
    end
end