function [count, unct] = compute_Moon_contour(lon,lat, params, Cnm, Snm)
    n_max = params(4);
    R_M = params(2)*1E3;    % [m]
    GM_moon = params(1);    % [m^3/s^2]
    normalized = params(3);
    h = 50E3;               % [m]

    ACAF_ACI = eye(3); 
    sigma = 1E-12;
    
    count = nan(length(lon(:, 1)), length(lon(1, :)));
    unct = nan(length(lon(:, 1)), length(lon(1, :)), 4);
    for k = 1:length(lon(1, :))
        for j = 1:length(lat(:, 1))
            % compute ECEF coordiantes
            phi = deg2rad(lat(j, k));   % [rad]
            lam = deg2rad(lon(j, k));   % [rad]

            x = (R_M).*cos(phi).*cos(lam);
            y = (R_M).*cos(phi).*sin(lam);
            z = (R_M).*sin(phi);
            r = [x;y;z];

            [~, ~, T_ECEF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                         r, R_M, GM_moon, normalized);
% %             [~, ~, T_ECEF0] = potentialGradient_nm(Cnm, Snm, 0, ...
% %                          r, R_M, GM_moon, normalized);
            ENU_ECEF = ecef2enu(r);
            T_ENU    = ENU_ECEF * T_ECEF * ENU_ECEF';

% %             DT_ECEF = T_ECEF - T_ECEF0;

            Y       = [T_ENU(1,1);T_ENU(1,2);T_ENU(2,2);T_ENU(3,3)]; 
            
            count(j, k) = norm(Y, 'fro');

            x = (R_M+h).*cos(phi).*cos(lam);
            y = (R_M+h).*cos(phi).*sin(lam);
            z = (R_M+h).*sin(phi);
            r = [x;y;z];
            
            ENU_ECEF = ecef2enu(r);
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, ...
                R_M, GM_moon, r, ACAF_ACI, ENU_ECEF');
            H = [Hpos(1:2, :); Hpos(5, :); Hpos(9, :)];   % Att - Less sensitive components
            % H = [Hpos(1:3, :); Hpos(5:6, :); Hpos(9, :)]; % Full Tensor
            % H = [Hpos(1, :); Hpos(5, :); Hpos(9, :)];     % Att and ang acc - Less sensitive components

            IF    = (1/sigma)^2 * (H' * H);
            P     = inv(IF);
            sg    = sqrt(diag(P));
            unct(j, k, 4) = sqrt(sg(1)^2 + sg(2)^2 + sg(3)^2);  % [m]
            unct(j, k, 2) = sg(2);  % [m]
            unct(j, k, 3) = sg(3);  % [m]
            unct(j, k, 1) = sg(1);  % [m]
        end
    end
end

