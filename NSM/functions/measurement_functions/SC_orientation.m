function [theta, M_ext] = SC_orientation(t, state_t, type, ACI_GRF, ECEF_ACI, Iner)
    % Description: compute S/C orientation over time given the time vector
    % and its position in the inertial frame. Theta is a 9 x Nt vector
    % which describes the attitude deviation from the inertial frame over
    % time.

    Nt  =length(t); 
    theta = ones(9, Nt) * NaN; % [theta; thetaDot; thetaDdot]
    R  = ones(3, 3, Nt) * NaN;
    if(type == "inertial")
        for j = 1:Nt
            R(: ,:, j) = eye(3,3);
        end
    elseif(type == "RTN")
        for j = 1:Nt
            r = state_t(j, 1:3)';    % [m]   ACI frame
            v = state_t(j, 4:6)';    % [m/s] ACI frame

            R_hat = r / vecnorm(r);
            N_hat = cross(r, v);
            N_hat = N_hat / vecnorm(N_hat);
            T_hat = cross(N_hat, R_hat);
        
            R(:, :, j) = [R_hat'; T_hat'; N_hat'];
        end
    elseif(type == "GRF")
        for j = 1:Nt
            maxPos = 3 * j; minPos = maxPos - 2;
            R(:, :, j) = ACI_GRF(minPos:maxPos, :); % [GRF_ACI matrix]
        end
    elseif(type == "ENU")
     for j = 1:Nt
        r = state_t(j, 1:3)';    % [m]   ACI frame

        maxVal = 3 *j; minVal = maxVal -2;
        ECEF_ACI_R = ECEF_ACI(minVal:maxVal, :);
        r_ECEF = ECEF_ACI_R * r;
        x = r_ECEF(1); y = r_ECEF(2); z = r_ECEF(3);

        % Longitude
        lambda = atan2(y, x);
        
        % Geocentric latitude
        phi = atan2(z, sqrt(x^2 + y^2));
        
        % Rotation matrix: ENU <- ECEF
        ENU_ECEF = [ -sin(lambda),              cos(lambda),               0;
                  -sin(phi)*cos(lambda),    -sin(phi)*sin(lambda),     cos(phi);
                   cos(phi)*cos(lambda),     cos(phi)*sin(lambda),     sin(phi) ];
        ENU_ACI  = ENU_ECEF * ECEF_ACI_R;
    
        R(:, :, j) = ENU_ACI;
    end

    else
        warning('Please, select a correct Instrument Frame (IF)')
    end

    % compute angular velocity (body frame)
    w = ones(3, Nt) * NaN; wd = w;
    for i = 2:Nt-1
        dt = t(i+1) - t(i-1);
        R1 = R(:,:,i-1);
        R2 = R(:,:,i+1);
        R0 = R(:,:,i);
        R_dot = (R2 - R1)./dt;
        omega_skew = -R_dot * R0';

        % Extract angular velocity from skew-symmetric matrix
        w(:,i) = [omega_skew(3,2); ...
            omega_skew(1,3); ...
            omega_skew(2,1)];
    end

    % Compute angular acceleration (body frame)
    for i = 2:Nt-1
        dt = t(i+1) - t(i-1);
        wd(:,i) = (w(:,i+1) - w(:,i-1)) / dt;
    end

    % compute external torque (body frame)
    M_ext  = ones(3, Nt) * NaN;
    for i = 1:Nt
        M_ext(:, i) = Iner*wd(:, i) + cross(w(:, i), Iner*w(:, i));
    end
    
    % Compute Euler angles and Euler angle rates 
    for j = 1:Nt
        BN = R(:, :, j);

        % compute yaw, pith, roll
        yaw   = atan2(BN(1,2), BN(1,1)); % psi
        pitch = -asin(BN(1,3));          % theta
        roll  = atan2(BN(2,3), BN(3,3)); % phi
        
        tol = 1E-4;
        if(abs(BN(1,3)) >= 1 - tol)
            roll = 0;
            yaw  = atan2(-BN(2,1), BN(2,2));
        end
    
        % store attitude angle
        theta(1:3, j) = [yaw;pitch;roll];
    
        % compute and store angle rate
        th = theta(2, j); phi = theta(3, j);
        T = (1/cos(th)).*[0, sin(phi), cos(phi);...
            0, cos(phi)*cos(th), -sin(phi)*cos(th);...
            cos(th), sin(phi)*sin(th), cos(phi)*sin(th)];
    
        theta(4:6, j) = T * w(:, j);
    end
end