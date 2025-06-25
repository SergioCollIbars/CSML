function [theta, M_ext] = SC_orientation(t, state_t, type, ACI_GRF, Iner)
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
% %         acc = ones(3, Nt) * NaN;
% %         v   = state_t(:, 4:6)';    % [m/s] ACI frame
% %         for j = 2:Nt-1
% %             dt = t(j+1) - t(j-1);
% %             acc(:, j) = (v(:, j+1) - v(:, j-1))./dt;
% %         end
% %         for j = 1:Nt
% %             r = state_t(j, 1:3)';    % [m]   ACI frame
% %             v = state_t(j, 4:6)';    % [m/s] ACI frame
% %             a = acc(:, j);           % [m/s^2] ACI frame
% % 
% %             R_hat = r / norm(r);
% %             N_hat = cross(r, v);
% %             N_hat = N_hat / norm(N_hat);
% %             T_hat = cross(N_hat, R_hat);
% %         
% %             BN = [R_hat'; T_hat'; N_hat'];
% % 
% %             % compute RTN frame derivatives
% %             h = cross(r, v);
% %             R_hat_dot = v/vecnorm(r) - r/(vecnorm(r)^3) * dot(r, v);
% %             N_hat_dot = cross(r, a)/vecnorm(h)  - h/(vecnorm(h)^2) * (dot(N_hat, cross(r, a)));
% %             T_hat_dot = cross(N_hat_dot, R_hat) + cross(N_hat, R_hat_dot);
% % 
% %             R_dot = [R_hat_dot'; T_hat_dot'; N_hat_dot'];
% %             wskew = - R_dot * BN';
% % 
% %             w(:, j) = [dot(N_hat,T_hat_dot);dot(R_hat, N_hat_dot);dot(T_hat, R_hat_dot)];
% % 
% %             % compute yaw, pith, roll
% %             yaw   = atan2(BN(1,2), BN(1,1)); % psi
% %             pitch = -asin(BN(1,3));          % theta
% %             roll  = atan2(BN(2,3), BN(3,3)); % phi
% %             
% %             tol = 1E-4;
% %             if(abs(BN(1,3)) >= 1 - tol)
% %                 roll = 0;
% %                 yaw  = atan2(-BN(2,1), BN(2,2));
% %             end
% % 
% %             % store attitude angle
% %             theta(1:3, j) = [yaw;pitch;roll];
% % 
% %             % compute and store angle rate
% %             th = theta(2, j); phi = theta(3, j);
% %             T = (1/cos(th)).*[0, sin(phi), cos(phi);...
% %                 0, cos(phi)*cos(th), -sin(phi)*cos(th);...
% %                 cos(th), sin(phi)*sin(th), cos(phi)*sin(th)];
% % 
% %             theta(4:6, j) = T * w(:, j);
% %         end 
% %         
% %         % get angular acceleration (finite differnces)
% %         wd = ones(3, Nt) * NaN;
% %         for j = 2:Nt-1
% %             dt = t(j+1) - t(j-1);
% %             wd(:, j) = (w(:, j+1) - w(:, j-1))./dt;
% %         end
% % 
% %         % compute and store angle acceleration
% %         for j = 1:Nt
% %             th = theta(2, j); phi = theta(3, j);
% %             A = [-sin(th), 0, 1;...
% %             sin(phi)*cos(th), cos(phi), 0;...
% %             cos(phi)*cos(th), -sin(phi), 0];
% % 
% %             thetaDot = theta(5, j); phiDot = theta(6, j);
% %              A_dot = [-cos(th)*thetaDot, 0, 0;...
% %             cos(phi)*cos(th)*phiDot-sin(phi)*sin(th)*thetaDot, -sin(phi)*phiDot, 0;...
% %             -sin(phi)*cos(th)*phiDot-cos(phi)*sin(th)*thetaDot, -cos(phi)*phiDot, 0];
% % 
% %             theta(7:9, j) = A\(wd(:, j) - A_dot * theta(4:6, j));
% %         end
    elseif(type == "GRF")
        for j = 1:Nt
            maxPos = 3 * j; minPos = maxPos - 2;
            R(:, :, j) = ACI_GRF(minPos:maxPos, :); % [GRF_ACI matrix]
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