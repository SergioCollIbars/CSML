function [X_final, P_final] = estimate_grav_field(n_solve, time, Nt, Y_true, ...
    rn, vn, planetParams, poleParams, Cnm, Snm, P0t, R, sigmaPos, sigmaBias)
    % extract asteroid parameteres
    GM = planetParams(1); Re = planetParams(2); n_max = n_solve;
    normalized = planetParams(4); 
    t = time;

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_solve); 
    [Nc_org, ~, ~] = count_num_coeff(planetParams(3)); 
    Pcc = P0t(2:Nc, 2:Nc);
    Pcs = P0t(2:Nc, Nc_org+1:Nc_org+Ns);
    Psc = P0t(Nc_org+1:Nc_org+Ns,2:Nc);
    Pss = P0t(Nc_org+1:Nc_org+Ns, Nc_org+1:Nc_org+Ns);

    P0 = [Pcc, Pcs;Psc, Pss];

    Pc  = diag([ones(6, 1)*(sigmaPos^4);sigmaBias'.^2]);
    Pxc = zeros(Ncs-1, 12);

    planetParams(3) = n_solve;

    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 10;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    while (count < iterMax) && (err < 1E3)
        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;

        [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)]./1E-12;
        
             % compute position and attitude partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos   = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;

            Hap    = compute_posPartials_2ndOrder(GM, ...
                rn_ACI(1), rn_ACI(2), rn_ACI(3));
            H2pos  = Hap./1E-12;

            H2pos_z  = -H2pos(1, :) - H2pos(4, :);
            H2pos_tot = [H2pos;H2pos_z];

            Hn = Hpos;

            % Compute NS
            c = null(Hn');
            if(isempty(c))
                [~,v,D] = svd(Hn');
                nv = length(diag(v));
                C = D(:, nv-1);
            else
                C = c;
            end
            
            % visibility matrix
            H = C' * Hc;
            Ha = C' * [H2pos_tot, eye(6)];
            r = C' * R * C;

            % prefit residuals
            dY = Y_true(:, j) - [Y_ACI(1:3); Y_ACI(5:6);Y_ACI(9)]./1E-12;
            dy = C' * dY;

            % solve normal equations
            ax = H' * (r\H);
            nx = H' * (r\dy);

            mxc = (H' * (r\Ha));
            mcc = (Ha'* (r\Ha));
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;

            Mxc = Mxc + mxc;
            Mcc = Mcc + mcc;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L(1:Ncs-1);
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  NSM update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    Px = inv(Ax_L);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pc*Sxc';

    P_final = Pxx(1:Ncs-1, 1:Ncs-1);
end

