function [R_new, Y_cal] = compute_meas_weigth(time,state_nom, ...
    Cnm, Snm, planetParams, poleParams, R0, Px,Y_cal)
    Nt = length(time);
    rn = state_nom(1:3, :);
    vn = state_nom(4:6, :);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    GM = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
    normalized = planetParams(4); 
    
% %     pref = nan(3, Nt);
    r     = sqrt(diag(R0(1, 1)));
    R_new = nan(6, 6, Nt);
    for j = 1:Nt
        rn_ACI = rn(:, j);

        p = reshape(Px(j, :), [12, 12]);
        sigmaP = sqrt(diag(p));
        dr = 5.*sigmaP(1:3);
        db = 5.*sigmaP(7:12);

        % ACAF to ACI rotation matrix
        Wt = W0 + W * time(j);
        ACAF_ACI  = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        B_ACI     = eye(3,3);
        ACAF_BODY = ACAF_ACI * B_ACI';

        % computed meas & select measurements
        [Y_ACI_1, ~] = gradiometer_meas(time(j) ,planetParams, ACAF_ACI,...
            [rn(:, j)', vn(:, j)'], zeros(9,1), Cnm, Snm);

        [Y_ACI_2, ~] = gradiometer_meas(time(j) ,planetParams, ACAF_ACI,...
            [(rn(:, j)+dr)', vn(:, j)'], zeros(9,1), Cnm, Snm);

        % first-order position partials
        [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
        Hpos   = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;

% %         C = null(Hpos');

% %         dY_bias = eye(6) * db;
        
        dY_pos = ([Y_ACI_2(1:3); Y_ACI_2(5:6);Y_ACI_2(9)] - ...
                  [Y_ACI_1(1:3); Y_ACI_1(5:6);Y_ACI_1(9)])./1E-12;
        dY_pos = dY_pos - (Hpos * dr);
        
% %         pref(:, j) = C' * (dY_pos + dY_bias);
        pref = abs((dY_pos));

        stats = pref;
        r_new = zeros(6,6);
        for m = 1:6
            r_new(m,m) = sqrt(r^2 + stats(m)^2);
        end
        R_new(:, :, j) = r_new.^2;
    end

    % compute std of the residuals
% %     stats = max(max(abs(pref)'));
% %     r     = sqrt(diag(R0(1, 1)));

% %     if(stats > r)
% %         noise = normrnd(0, stats, [6, Nt]);
% %         r_new = sqrt(r^2 + stats^2);
% %         R_new = (1).*diag(ones(6, 1).*r_new).^2;
% %     else
% %         noise = zeros(6, Nt);
% %         R_new  = R0;
% %     end

% %     noise = normrnd(0, stats, [6, Nt]);
% %     r_new = sqrt(r^2 + stats^2);
% %     R_new = (1).*diag(ones(6, 1).*r_new).^2;
% % 
% %     Y_cal = Y_cal + noise;
end

