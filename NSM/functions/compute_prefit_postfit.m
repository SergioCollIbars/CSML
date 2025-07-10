function [prefit,postfit] = compute_prefit_postfit(data, Cnm, Snm,  dx, ...
    RotPlanet, B_ACI_mat, attitude, rn, vn, angVel, angAcc, t, planetParams)
    Nt  = length(data(1, :)); 
    prefit = zeros(9, Nt); postfit = prefit;

    for j = 1:Nt
        % Planet orientation
        maxPos = 3*j; minPos = maxPos - 2;
        ACAF_ACI = RotPlanet(minPos:maxPos, :);
    
        % from ACI to Nominal body frame
        B_ACI = B_ACI_mat(minPos:maxPos, :);

        [Y_ACI, Hc_ACI, ~] = gradiometer_meas(t(j) ,planetParams, ACAF_ACI, [rn(:, j)', vn(:, j)'], ...
                    zeros(9, Nt), Cnm, Snm, B_ACI');
        Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);

        [Yc] = add_angularComponents(Y_ACI, attitude(:, j), zeros(3, Nt), angVel(:, j),...
            angAcc(:, j));
        
        % extract angular acceleration effects
        D = [data(1, j),data(2, j),data(3, j);data(4, j),data(5, j),data(6, j);...
            data(7, j),data(8, j),data(9, j)];
        OMEGA_c = 0.5 * (D - D');
        omegac = [OMEGA_c(1,1);OMEGA_c(1,2);OMEGA_c(1,3);OMEGA_c(2,1);...
            OMEGA_c(2,2);OMEGA_c(2,3);OMEGA_c(3,1);OMEGA_c(3,2);OMEGA_c(3,3)];
        omegac = zeros(9, 1);

        % prefit 
        Y = data(:, j) - omegac;
        prefit(:, j)  = Y - Yc;
        postfit(:, j) = prefit(:, j) - Hc_BODY * dx; 
    end
end

