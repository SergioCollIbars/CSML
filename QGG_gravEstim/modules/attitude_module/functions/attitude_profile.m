function [phi, phiDot, phiDdot, theta,thetaDot, thetaDdot, psi, psiDot, psiDdot] = attitude_profile(C1, C2, C3, C4, C5, axis, t)
      %%                   ATTITUDE PROFILE FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This function is in charge of compute the spacecraft 
    %    profile given a attitude profile
    %
    %   Input:
    %       C1: first constant term
    %       C2: second constant term
    %       C3: third constant term
    %       C4: fourth constant term
    %       C5: fith constant term
    %       axis: array with the axis affected
    %       t: time vector
    %       shift: value shift
    %   
    %   Output: 
    %       phi: roll angle [rad]
    %       phiDot: roll angular velocity [rad / s]
    %       phiDdot: roll angular acc [rad/ s^2]
    %       theta: ptich angle [rad]
    %       thetaDot: ptich angular velocity [rad / s]
    %       thetaDdot: ptich angular acc [rad/ s^2]
    %       psi: yaw angle [rad]
    %       psiDot: psi angular velocity [rad / s]
    %       psiDdot: psi angular acc [rad/ s^2]
    % --------------------------------------------------------------------%
    % get axis affected
    axisN = length(axis);
    axis1 = 0;
    axis2 = 0;
    axis3 = 0;

    for i = 1: axisN
        if(axis(i) == 1)
            axis1 = 1;
        elseif(axis(i) == 2)
            axis2 = 1;
        elseif(axis(i) == 3)
            axis3 = 1;
        end
    end

    % Euler angles matrix definition
    phi = zeros(1, length(t));          % R[AXIS 1]
    theta = zeros(1, length(t));        % R[AXIS 2]
    psi = zeros(1, length(t));          % R[AXIS 3]
    
    phiDot = zeros(1, length(t));       % R[AXIS 1]
    thetaDot = zeros(1, length(t));     % R[AXIS 2]
    psiDot = zeros(1, length(t));       % R[AXIS 3]
    
    phiDdot = zeros(1, length(t));      % R[AXIS 1]
    thetaDdot = zeros(1, length(t));    % R[AXIS 2]
    psiDdot = zeros(1, length(t));      % R[AXIS 3]
    
    % compute attitude profile
    for k =1:length(t)
        if(axis1 == 1)
            phi(k)      = C1 + C2 * t(k) + C3 * sin(C4 * t(k) + C5);
            phiDot(k)   = C2 + C3 * C4 * cos(C4 * t(k) + C5);
            phiDdot(k)  = - C3 * C4 * C4 * sin(C4 * t(k) + C5);
        end
        if(axis2 == 1)
            theta(k)    = C1 + C2 * t(k) + C3 * sin(C4 * t(k) + C5);
            thetaDot(k) = C2 + C3 * C4 * cos(C4 * t(k) + C5);
            thetaDdot(k)= - C3 * C4 * C4 * sin(C4 * t(k) + C5);
        end
        if(axis3 == 1)
            psi(k)      = C1 + C2 * t(k) + C3 * sin(C4 * t(k) + C5);
            psiDot(k)   = C2 + C3 * C4 * cos(C4 * t(k) + C5);
            psiDdot(k)  = - C3 * C4 * C4 * sin(C4 * t(k) + C5);
        end
    end
        
end
