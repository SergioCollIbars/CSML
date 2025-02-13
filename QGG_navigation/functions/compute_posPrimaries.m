function [posE, posM, posS] = compute_posPrimaries(TIME, planetParams, system)

    posE = ones(3, length(TIME)) * NaN;
    posM = posE; posS = posE;

    for j = 1:length(TIME)
        % current time
        t = TIME(j);    % [-]

        if(system == "EPHEM")
            [pos_earth, pos_moon, pos_sun] = computePos_SPICE(t, planetParams);
        elseif(system == "CR3BP" || system == "FCR3BP")
            [pos_earth, pos_moon, pos_sun] = computePos_circular(t, planetParams);
        end
        
        % save values
        posE(:, j) = pos_earth;  % [-]
        posM(:, j) = pos_moon;   % [-]
        posS(:, j) = pos_sun;    % [-]
    end
end


function [pos_earth, pos_moon, pos_sun] = computePos_SPICE(t, planetParams)
    % Earth position
    target = 'EARTH';
    et = t./planetParams(3); % Convert UTC time to ephemeris time
    ref = 'J2000';
    abcorr = 'NONE';
    observer = '3';          % Set the observer to the Earth-Moon barycenter

    [pos1, ~] = cspice_spkpos(target, et, ref, abcorr, observer);
    pos1 = pos1./planetParams(2) * 1E3;
    pos_earth = pos1;   % [-]
    
    % Moon position
    target = 'MOON';
    [pos2, ~] = cspice_spkpos(target, et, ref, abcorr, observer);
    pos2 = pos2./planetParams(2) * 1E3;
    pos_moon = pos2;    % [-]

    % Sun position
    target = 'SUN';
    [pos3, ~] = cspice_spkpos(target, et, ref, abcorr, observer);
    pos3 = pos3./planetParams(2) * 1E3;
    pos_sun = pos3;    % [-]
end

function [pos_earth, pos_moon, pos_sun] = computePos_circular(t, planetParams)
    mu = planetParams(1);
    M  = t;
    ns = 2*pi/(365*86400);                  % [rad/s]
    AU = 1.496e+11 / planetParams(2);       % [-]

    % Earth position
    pos_earth = -[mu*cos(M);mu*sin(M);0];

    % Moon position
    pos_moon  = -[(mu-1)*cos(M);(mu-1)*sin(M); 0];

    % Sun position
    scale = ns/planetParams(3);
    pos_sun   = [AU*cos(M*scale+pi/4); AU*sin(M*scale+pi/4); 0];  
end
