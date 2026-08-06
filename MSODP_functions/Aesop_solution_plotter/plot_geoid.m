function [] = ...
    plot_geoid(Cnm, Snm, RE, gridResolution, tt, limVec)
    %PLOT_GEOID_AND_WATER_HEIGHT
    % Synthesize and plot:
    %   1. Geoid-height error [m]
    %   2. Equivalent water height [cm]
    %
    % MATRIX CONVENTION
    %   Cnm(n+1,m+1) = Delta C_nm
    %   Snm(n+1,m+1) = Delta S_nm
    %
    % Entries for m > n may be NaN.
    %
    % INPUTS
    %   Cnm, Snm       Spherical-harmonic coefficient-error matrices
    %                  Size: (Nmax+1)-by-(Nmax+1)
    %
    %   RE             Earth reference radius [m]
    %                  Default: 6378136.3 m
    %
    %   gridResolution Latitude/longitude grid resolution [deg]
    %                  Default: 1 deg
    %
    %   tt             Graphic title (string)
    %
    %   limVec         Color limits vector [min max]
    %
    % OUTPUTS
    %   N/A
    %
    % NOTES
    %   - Coefficients are assumed fully normalized.
    %   - Degrees 0 and 1 are ignored.
    %   - Snm(n+1,1), corresponding to S_n0, is set to zero because
    %     zonal sine coefficients do not contribute.

    %% Default inputs

    if nargin < 3 || isempty(RE)
        RE = 6378136.3;           % [m]
    end

    if nargin < 4 || isempty(gridResolution)
        gridResolution = 1;       % [deg]
    end

    Nmax = size(Cnm,1) - 1;

    %% Input checks

    if ~isnumeric(Cnm) || ~isnumeric(Snm)
        error('Cnm and Snm must be numeric matrices.');
    end

    if ~isequal(size(Cnm),size(Snm))
        error('Cnm and Snm must have the same dimensions.');
    end

    if size(Cnm,1) ~= size(Cnm,2)
        error('Cnm and Snm must be square matrices.');
    end

    %% Check coefficients in the valid triangular region

    for n = 2:Nmax

        for m = 0:n

            C = Cnm(n+1,m+1);

            if ~isfinite(C)
                error(['Cnm(%d,%d), corresponding to C_%d,%d, ', ...
                       'is not finite.'], ...
                       n+1,m+1,n,m);
            end

            % S_n0 does not contribute and may be stored as NaN
            if m > 0

                S = Snm(n+1,m+1);

                if ~isfinite(S)
                    error(['Snm(%d,%d), corresponding to S_%d,%d, ', ...
                           'is not finite.'], ...
                           n+1,m+1,n,m);
                end

            end

        end

    end

    %% Latitude-longitude grid

    latDeg = -90:gridResolution:90;
    lonDeg = -180:gridResolution:180;

    latRad = deg2rad(latDeg);
    lonRad = deg2rad(lonDeg);

    nLat = length(latDeg);
    nLon = length(lonDeg);

    geoidSum = zeros(nLat,nLon);

    % Associated Legendre argument
    x = sin(latRad);
    
    %% Spherical-harmonic synthesis

    for n = 2:Nmax

        % Pbar(m+1,ilat) = Pbar_nm(sin(latitude))
        Pbar = fully_normalized_legendre(n,x);

        for m = 0:n

            C = Cnm(n+1,m+1);

            if m == 0
                % S_n0 is always zero
                S = 0;
            else
                S = Snm(n+1,m+1);
            end

            if C == 0 && S == 0
                continue
            end

            longitudeTerm = ...
                C*cos(m*lonRad) + ...
                S*sin(m*lonRad);

            degreeOrderTerm = ...
                Pbar(m+1,:).' * longitudeTerm;

            % Dimensionless geoid SH sum
            geoidSum = geoidSum + degreeOrderTerm;

        end

    end

    %% Convert to physical units

    % Geoid height [m]
    geoidError = RE*geoidSum;

    %% Plot maps

    figure;

    % Geoid-height error
    nexttile;

    imagesc(lonDeg,latDeg,geoidError);

    set(gca,'YDir','normal');
    axis equal tight;

    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');

    title(sprintf( ...
        'Geoid-height error up to degree %d',Nmax));

    cb1 = colorbar;
    cb1.Label.String = 'Geoid-height error';

    colormap(turbo);

    hold on;
    plot_coastlines();
    hold off;

    title(tt);
    clim(limVec);

end


function Pbar = fully_normalized_legendre(n,x)
%FULLY_NORMALIZED_LEGENDRE
% Compute fully normalized associated Legendre functions:
%
%   Pbar_nm =
%   sqrt[(2-delta_m0)(2n+1)(n-m)!/(n+m)!] P_nm
%
% MATLAB's legendre function includes the Condon-Shortley phase
% (-1)^m, which is removed here.

    P = legendre(n,x);

    Pbar = zeros(size(P));

    for m = 0:n

        factorialRatio = exp( ...
            gammaln(n-m+1) - gammaln(n+m+1));

        if m == 0
            orderFactor = 1;
        else
            orderFactor = 2;
        end

        normalization = sqrt( ...
            orderFactor*(2*n+1)*factorialRatio);

        % Remove Condon-Shortley phase
        Pbar(m+1,:) = ...
            (-1)^m*normalization*P(m+1,:);

    end

end


function plot_coastlines()
%PLOT_COASTLINES Plot MATLAB coastline data when available.

    try

        coast = load('coastlines');

        plot( ...
            coast.coastlon, ...
            coast.coastlat, ...
            'k', ...
            'LineWidth',0.7);

    catch

        % Continue without coastlines if the file is unavailable

    end

end