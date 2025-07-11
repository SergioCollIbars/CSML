function [angVel, angAcc] = compute_angularVals_v2(theta, thetaDot, thetaDdot)
    % Description: given the Euler angles in a (3-2-1) rotation and the
    % Euler angles rate (velocity and acceleration) compute the angular
    % velocity and acceleration in the body frame. 
    % Compute the partials of the angular values with respect to the Euler
    % angles.
    Nt = length(theta(1, :));
    angAcc = zeros(3, Nt);
    dA_dt    = thetaDot;
    attitude = theta;
    vals = ones(3, Nt) * NaN;
    for j = 1:Nt                            % attitude rate of change (vel)
        theta = attitude(2, j); phi = attitude(3, j);
        A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
        
        % angular velocity
        vals(:, j) = A * dA_dt(: ,j);
    end
    angVel = vals;

    for j = 1:Nt
        theta = attitude(2, j); phi = attitude(3, j);
        psiDot = dA_dt(1, j); thetaDot = dA_dt(2, j); phiDot = dA_dt(3, j);
        A = [-sin(theta), 0, 1;...
            sin(phi)*cos(theta), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0];
        A_dot = [-cos(theta)*thetaDot, 0, 0;...
            cos(phi)*cos(theta)*phiDot-sin(phi)*sin(theta)*thetaDot, -sin(phi)*phiDot, 0;...
            -sin(phi)*cos(theta)*phiDot-cos(phi)*sin(theta)*thetaDot, -cos(phi)*phiDot, 0];

        angAcc(:, j) = A_dot * dA_dt(:, j) + A * thetaDdot(:, j) ;
    end
end