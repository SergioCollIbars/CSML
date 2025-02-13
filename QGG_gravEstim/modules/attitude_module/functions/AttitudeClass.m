classdef AttitudeClass
    %%                      ATTITUDE CLASS
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 26/10/2022
    %
    %   Description: Class to store all the information related with the
    %   attitude module
    %
    % --------------------------------------------------------------------%
    properties
        % profile configuration
        AttitudeData

        % First rotation angles
        theta1
        theta1Dot
        theta1Ddot

        % second rotation angles
        theta2
        theta2Dot
        theta2Ddot

        % third rotation nagles
        theta3
        theta3Dot
        theta3Ddot

        % angular velocity
        omega
        %angular acceleration
        Omega

        % time vector
        t
    end

%     methods
%         function obj = untitled2(inputArg1,inputArg2)
%             %UNTITLED2 Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
% 
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
%     end
end