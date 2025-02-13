function [AttitudeObj] = angular_vel_acc(AttitudeObj, t)
      %%            ANGULAR VELOCITY & ACCELERATION FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This function is in charge of compute the angular
    %   velocities and accelerations given an attitude profile and an Euler
    %   rotation order.
    %
    %   Input:
    %       AttitudeObj: object to staore attitude varaibles
    %       t: time vector
    %
    %   Output: 
    %       AttitudeObj: object to staore attitude varaibles
    %    
    % --------------------------------------------------------------------%

    % get angle first roation
    theta1 = AttitudeObj.theta1;
    theta1Dot = AttitudeObj.theta1Dot;
    theta1Ddot = AttitudeObj.theta1Ddot;
    
    % get angle second rotation
    theta2 = AttitudeObj.theta2;
    theta2Dot = AttitudeObj.theta2Dot;
    theta2Ddot = AttitudeObj.theta2Ddot;

    % get angle third rotation
    theta3 = AttitudeObj.theta3;
    theta3Dot = AttitudeObj.theta3Dot;
    theta3Ddot = AttitudeObj.theta3Ddot;

    % matrix definition
    omega = zeros(3, length(t));
    Omega = zeros(3, length(t));

    for k = 1:length(t)
        % angular velocities
        omega(1, k) = theta1Dot(k) * sin(theta3(k)) * sin(theta2(k)) + ...
            cos(theta3(k)) * theta2Dot(k);
        omega(2, k) = cos(theta3(k)) * sin(theta2(k)) * theta1Dot(k) - ...
            sin(theta3(k)) * theta2Dot(k);
        omega(3, k) = theta1Dot(k) * cos(theta2(k)) + theta3Dot(k);

        % angular acceleration
        Omega(1, k) = theta1Ddot(k) * (sin(theta3(k)) * sin(theta2(k))) + ...
            theta1Dot(k) * (cos(theta3(k)) * sin(theta2(k)) * theta3Dot(k) ...
            + sin(theta3(k)) * cos(theta2(k)) * theta2Dot(k)) ...
            + theta2Ddot(k) * cos(theta3(k)) - theta2Dot(k) * ...
            sin(theta3(k)) * theta3Dot(k);
        Omega(2, k) = theta1Ddot(k) * cos(theta3(k)) * sin(theta2(k)) + ...
            theta1Dot(k) * (-sin(theta3(k)) * theta3Dot(k) * sin(theta2(k)) ...
            + cos(theta3(k)) * cos(theta2(k)) * theta2Dot(k))...
            - theta2Ddot(k) * sin(theta3(k)) - theta2Dot(k) * cos(theta3(k)) ...
            * theta3Dot(k);
        Omega(3, k) = theta1Ddot(k) * cos(theta2(k)) - theta1Dot(k) * ...
            sin(theta2(k)) * theta2Dot(k) + theta3Ddot(k);

    end

    % store variables into object
    AttitudeObj.omega = omega;
    AttitudeObj.Omega = Omega;

end