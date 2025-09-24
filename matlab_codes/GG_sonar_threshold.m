clear;
clc;
close all;

format long g;

%%              GRADIOMETER SONAR

% constants
G = 6.67E-11;                       % [m^3 Kg^-1 s^-2]

% mass & velocity range 
N = 1000;
M_range = linspace(1, 5E3, N);    % [Kg]
V_range = linspace(10, 50E3, N);  % [m/s]

% Gradiometer measurement noise STD & Sampling interval
At = 5;                             % [sec]
sigma_m = sqrt(2) * 1E-15;          % [s^-2]
T_dot = sigma_m / At;               % [s^-3]
T = sigma_m / sqrt(2);


% Maximum radius
r_max = ones(length(V_range), length(M_range)) * NaN;
r_max_2 = ones(1, length(M_range)) * NaN;
for i = 1:N
    V = V_range(i);
    for j = 1:N
        M = M_range(j);
        r_max(i, j) = (18 * G * M / T_dot * V * cos(pi/4))^(0.25);
    end
end

for k = 1:N
    M = M_range(k);
    r_max_2(k) = (6 * G * M / T)^(1/3);
end

% plot
plotVals(M_range, V_range, r_max, r_max_2)

figure()
plot(M_range, r_max_2, 'LineWidth', 2);
xlabel('Mass [Kg]')
ylabel('[m]')
title('Maximum range based on gradiometer signal')
grid on;

%% FUNCTIONS
function [] = plotVals(xv, yv, F, x)
% Inputs:
% F  : Ny x Nx   (Z values)
% xv : 1 x Nx    (x-coordinates for columns)
% yv : 1 x Ny    (y-coordinates for rows)
% x  : 1 x Nx    (target Z level per column)

figure
contourf(xv, yv./1E3, F, 40, 'LineStyle', 'none');  % filled only
hold on
colorbar
colormap('turbo')
clim([0 1500])     % set color scale limits
xlabel('Mass  [Kg]'); ylabel('Velocity [Km / s]'); 
title('Maximum range based on gradiometer signal rate')
set(gca,'XScale','log','YScale','log')      % log scale for both axes

Ny = size(F,1);
Nx = size(F,2);
y_curve = nan(1, Nx);

for j = 1:Nx
    zlev = x(j);                 % desired Z level in this column
    col  = F(:,j);               % Ny×1
    d    = col - zlev;

    % find sign changes (crossings) between consecutive samples
    sgn      = d(1:end-1).*d(2:end);
    k_cross  = find(sgn <= 0 & ~isnan(sgn));   % indices where a crossing occurs between k and k+1

    if isempty(k_cross)
        % no crossing found in this column -> leave NaN
        continue
    end

    % linearly interpolate each crossing to get a y-value
    ylist = zeros(numel(k_cross),1);
    for k = 1:numel(k_cross)
        k0 = k_cross(k);
        z0 = col(k0);
        z1 = col(k0+1);
        y0 = yv(k0);
        y1 = yv(k0+1);
        if z1 ~= z0
            t = (zlev - z0) / (z1 - z0);
            ylist(k) = y0 + t*(y1 - y0);
        else
            ylist(k) = (y0 + y1)/2;  % flat segment exactly at level
        end
    end

    % If multiple crossings, choose the one nearest the column maximum
    [~, imax] = max(col);
    % map each crossing to the nearest row index to compare proximity to imax
    [~, k_near] = min(abs(k_cross - imax));
    y_curve(j) = ylist(k_near);
end

% Overlay the curve
loglog(xv, y_curve./1E3,'LineWidth', 2, 'LineStyle','--', 'Color','k')
grid on;

% (Optional) log axes if needed:
% set(gca, 'XScale','log', 'YScale','log')

end