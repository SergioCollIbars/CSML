clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%% PLOT TIME VARIABLE GRAV. FIELD GEOID
% Description: using the GOCO2025 gravity model downloaded from:
%   https://icgem.gfz.de/home. Compute the GEOID heigh in mm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'GOCO2025s.gfc';

Lmax = 120;

% Evaluation date
% Change this to the epoch you want.
date_eval = datetime(2008,8,31);

% Reference epoch from GOCO2025s header
date_ref = datetime(2015,1,1);

% Annual period [days]
T = 365.25;

% Time difference in years
dt_years = days(date_eval - date_ref)/T;

% Earth radius from GOCO2025s header [m]
R = 6.3781363000e6;

% Allocate coefficient arrays
C0      = zeros(Lmax+1,Lmax+1);
S0      = zeros(Lmax+1,Lmax+1);

Ctrnd   = zeros(Lmax+1,Lmax+1);
Strnd   = zeros(Lmax+1,Lmax+1);

Cacos1  = zeros(Lmax+1,Lmax+1);
Sacos1  = zeros(Lmax+1,Lmax+1);
Casin1  = zeros(Lmax+1,Lmax+1);
Sasin1  = zeros(Lmax+1,Lmax+1);

Cacos05 = zeros(Lmax+1,Lmax+1);
Sacos05 = zeros(Lmax+1,Lmax+1);
Casin05 = zeros(Lmax+1,Lmax+1);
Sasin05 = zeros(Lmax+1,Lmax+1);

fid = fopen(filename,'r');

if fid < 0
    error('Could not open file: %s', filename);
end

while ~feof(fid)

    line = strtrim(fgetl(fid));

    if isempty(line)
        continue
    end

    parts = strsplit(line);
    key = parts{1};

    if ~ismember(key, {'gfc','gfct','trnd','acos','asin'})
        continue
    end

    l = str2double(parts{2});
    m = str2double(parts{3});

    if l > Lmax
        continue
    end

    C = str2double(parts{4});
    S = str2double(parts{5});

    ii = l + 1;
    jj = m + 1;

    switch key

        case {'gfc','gfct'}
            C0(ii,jj) = C;
            S0(ii,jj) = S;

        case 'trnd'
            Ctrnd(ii,jj) = C;
            Strnd(ii,jj) = S;

        case 'acos'
            period = str2double(parts{8});

            if abs(period - 1.0) < 1e-12
                Cacos1(ii,jj) = C;
                Sacos1(ii,jj) = S;
            elseif abs(period - 0.5) < 1e-12
                Cacos05(ii,jj) = C;
                Sacos05(ii,jj) = S;
            end

        case 'asin'
            period = str2double(parts{8});

            if abs(period - 1.0) < 1e-12
                Casin1(ii,jj) = C;
                Sasin1(ii,jj) = S;
            elseif abs(period - 0.5) < 1e-12
                Casin05(ii,jj) = C;
                Sasin05(ii,jj) = S;
            end
    end
end

fclose(fid);

% Time-variable coefficient anomaly relative to t0
dC = Ctrnd  * dt_years ...
   + Cacos1 * cos(2*pi*dt_years) ...
   + Casin1 * sin(2*pi*dt_years) ...
   + Cacos05 * cos(4*pi*dt_years) ...
   + Casin05 * sin(4*pi*dt_years);

dS = Strnd  * dt_years ...
   + Sacos1 * cos(2*pi*dt_years) ...
   + Sasin1 * sin(2*pi*dt_years) ...
   + Sacos05 * cos(4*pi*dt_years) ...
   + Sasin05 * sin(4*pi*dt_years);

% Degree geoid-height amplitude
% N_l = R * sqrt(sum_m dC_lm^2 + dS_lm^2)
degree = (0:Lmax)';
N_l_m  = zeros(Lmax+1,1);

for l = 0:Lmax
    ii = l + 1;
    jj = 1:(l+1);

    N_l_m(ii) = R * sqrt( ...
        sum(dC(ii,jj).^2 + dS(ii,jj).^2) );
end

% Convert to mm
N_l_mm = 1000 * N_l_m;

% Main plot
figure;
semilogy(degree, N_l_mm, 'LineWidth', 1.5, 'Color', 'k');
grid on;
xlabel('Spherical harmonic degree');
ylabel('Geoid height [mm]');
title(['GOCO2025s time-variable geoid signal on ', datestr(date_eval)]);
xlim([1 Lmax]);   % start at 1 because spatial resolution is undefined for n = 0

ax1 = gca;

% Choose degree ticks for the main axis
deg_ticks = [20 40 60 80 100 120];
deg_ticks = deg_ticks(deg_ticks <= Lmax);

ax1.XTick = deg_ticks;
ax1.XLim  = [1 Lmax];

% Create a second x-axis below the main one
ax2 = axes( ...
    'Position', ax1.Position, ...
    'XAxisLocation', 'bottom', ...
    'YAxisLocation', 'right', ...
    'Color', 'none', ...
    'YColor', 'none', ...
    'XColor', 'k', ...
    'XLim', ax1.XLim, ...
    'YLim', ax1.YLim);

% Move the second axis slightly downward
ax2.Position(2) = ax1.Position(2) - 0.08;

% Use same x tick locations, but label them as spatial resolution
ax2.XTick = deg_ticks;

spatial_res_km = 40000 ./ (2 * deg_ticks);

ax2.XTickLabel = compose('%.0f', spatial_res_km);
xlabel(ax2, 'Approximate spatial resolution [km]');

% Keep the original axis on top
uistack(ax1, 'top');

% Make sure both axes stay linked
linkaxes([ax1 ax2], 'x');