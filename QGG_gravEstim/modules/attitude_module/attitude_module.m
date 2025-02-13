function [AttitudeObj] = attitude_module(AttitudeObj, t, save)
    %%                      ATTITUDE MODULE
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: This module is in charge of compute the angular
    %   velocities and accelerations given an attitude profile and an Euler
    %   rotation order.
    %
    %   Input:
    %       AttitudeObj:  attitude class object
    %       t: time vector
    %       save: save variables option boolean
    %
    %   Output:
    %       AttitudeObj:  attitude class object
    % --------------------------------------------------------------------%
    
    % add path
    addpath('modules/orbital_module/functions/')

    disp('  Attitude module')
    
    tic
    % first rotation
    C1  = deg2rad(AttitudeObj.AttitudeData(1));
    C2  = deg2rad(AttitudeObj.AttitudeData(2));
    C3  = deg2rad(AttitudeObj.AttitudeData(3));
    C4  = deg2rad(AttitudeObj.AttitudeData(4));
    C5  = deg2rad(AttitudeObj.AttitudeData(5));

    [~, ~, ~, ~, ~, ~, AttitudeObj.theta1, AttitudeObj.theta1Dot, AttitudeObj.theta1Ddot] = ...
        attitude_profile(C1, C2, C3, C4, C5, [0, 0, 3], t);

    % second rotation
    C1  = deg2rad(AttitudeObj.AttitudeData(6));
    C2  = deg2rad(AttitudeObj.AttitudeData(7));
    C3  = deg2rad(AttitudeObj.AttitudeData(8));
    C4  = deg2rad(AttitudeObj.AttitudeData(9));
    C5  = deg2rad(AttitudeObj.AttitudeData(10));

    [AttitudeObj.theta2, AttitudeObj.theta2Dot, AttitudeObj.theta2Ddot, ~, ~, ~, ~, ~, ~] = ...
        attitude_profile(C1, C2, C3, C4, C5, [1, 0, 0], t);

    % third roatation
    C1  = deg2rad(AttitudeObj.AttitudeData(11));
    C2  = deg2rad(AttitudeObj.AttitudeData(12));
    C3  = deg2rad(AttitudeObj.AttitudeData(13));
    C4  = deg2rad(AttitudeObj.AttitudeData(14));
    C5  = deg2rad(AttitudeObj.AttitudeData(15));

    [~, ~, ~, ~, ~, ~, AttitudeObj.theta3, AttitudeObj.theta3Dot, AttitudeObj.theta3Ddot] = ...
        attitude_profile(C1, C2, C3, C4, C5, [0, 0, 3], t);

    % compute angular velocity
    AttitudeObj = angular_vel_acc(AttitudeObj, t);
    toc
    
    % save time
    AttitudeObj.t = t;
    
    % save values if true
    if(save == true)
        saveData_attitude(AttitudeObj, "attitudeData.txt")
    end

end