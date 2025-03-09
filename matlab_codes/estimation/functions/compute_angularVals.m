function [angVel, angAcc, H_angVel, H_angAcc] = compute_angularVals(attitude, dA_dt, ddA_ddt, t)
    % Description: given the Euler angles in a (3-2-1) rotation and the
    % Euler angles rate (velocity and acceleration) compute the angular
    % velocity and acceleration in the body frame. 
    % Compute the partials of the angular values with respect to the Euler
    % angles.
    Nt = length(t);
    angVel = zeros(3, Nt); angAcc = zeros(3, Nt);
    dt = t(2) - t(1);                       % WARNING: assumes equal time
    for j = 1:Nt                            % attitude rate of change (vel)
        theta = attitude(2, j); phi = attitude(3, j);
        A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
        
        % angular velocity
        angVel(:, j) = A * dA_dt(: ,j);
        
        % partials of angular velocity w.r.t Euler angles
        finish = 3*j;
        init = finish - 2;
        psiDot = dA_dt(1, j); thetaDot = dA_dt(2, j);
        B = [0, -cos(theta)*psiDot, 0;...
             0, -sin(phi)*sin(theta)*psiDot, cos(phi)*cos(theta)*psiDot-sin(phi)*thetaDot;...
             0, -cos(phi)*sin(theta)*psiDot, -sin(phi)*cos(theta)*psiDot-cos(phi)*thetaDot];
        H_angVel(init:finish, :) = B + (1/dt).*A; 
    end
    for j = 1:Nt
        theta = attitude(2, j); phi = attitude(3, j);
        psiDot = dA_dt(1, j); thetaDot = dA_dt(2, j); phiDot = dA_dt(3, j);
        psiDdot = ddA_ddt(1, j); thetaDdot = ddA_ddt(2, j); phiDdot = ddA_ddt(3, j);
        A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
        A_dot = [-cos(theta)*thetaDot, 0, 0;...
            cos(phi)*cos(theta)*phiDot-sin(phi)*sin(theta)*thetaDot, -sin(phi)*phiDot, 0;...
            -sin(phi)*cos(theta)*phiDot-cos(phi)*sin(theta)*thetaDot, -cos(phi)*phiDot, 0];

        % angular acceleration
        angAcc(:, j) = A_dot * dA_dt(:, j) + A * ddA_ddt(:, j);

        % partials of angular acceleration w.r.t Euler angles
        C12 = (sin(theta)*thetaDot-cos(theta)/dt)*psiDot;
        C22 = (-cos(phi)*sin(theta)*phiDot - sin(phi)*cos(theta)*thetaDot - sin(phi)*sin(theta)/dt)*psiDot;
        C31 = (sin(phi)*sin(theta)*phiDot -cos(phi)*cos(theta)*thetaDot - cos(phi)*sin(theta)/dt)*psiDot;
        C23 = (-sin(phi)*cos(theta)*phiDot + cos(phi)*cos(theta)/dt - cos(phi)*sin(theta)*thetaDot)*psiDot + ...
            (-cos(phi)*phiDot -sin(phi)/dt)*thetaDot;
        C33 = (-cos(phi)*cos(theta)*phiDot - sin(phi)*cos(theta)/dt + sin(phi)*sin(theta)*thetaDot)*psiDot + ...
            (sin(phi)*phiDot -cos(phi)/dt)*thetaDot;
        C = [0,C12,0;...
             0,C22,C23;...
             0,C31,C33];

        D = [0, -cos(theta)*psiDdot, 0;...
             0, -sin(phi)*sin(theta)*psiDdot, cos(phi)*cos(theta)*psiDdot-sin(phi)*thetaDdot;...
             0, -cos(phi)*sin(theta)*psiDdot, -sin(phi)*cos(theta)*psiDdot-cos(phi)*thetaDdot];
        finish = 3*j;
        init = finish - 2;
        H_angAcc(init:finish, :) = (1/dt).*A_dot + (1/dt^2).*A + C + D;
    end
end
