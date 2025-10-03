function [R] = rotationMatrix(theta1, theta2, theta3, axis)
    %%                     ROTATION MATRIX FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 31/10/2022
    %
    %   Description: This function computes the rotation matrix for a given
    %       Euler angles order. 
    %
    %   Input:
    %       theta1: first angle rotation
    %       theta2: second angle rotation
    %       theta3: third angle rotation
    %       axis: axis matrix with rotation order
    %
    %   Output:
    %       OrbitObj:  orbit class object
    % --------------------------------------------------------------------%

     % matrix definition
    R1 = zeros(3,3);
    R2 = zeros(3,3);
    R3 = zeros(3,3);
    
    % theta 1.
    if(axis(1) == 1)
        R1 = [1, 0, 0;
              0, cos(theta1), sin(theta1);
              0, -sin(theta1), cos(theta1)];
    elseif(axis(1) == 2)
        R1 = [cos(theta1), 0, -sin(theta1);
              0, 1, 0;
              sin(theta1), 0, cos(theta1)];
    elseif(axis(1) == 3)
        R1 = [cos(theta1), sin(theta1), 0;                       
              -sin(theta1), cos(theta1), 0;
              0, 0, 1];
    end

    % theta 2.
    if(axis(2) == 1)
        R2 = [1, 0, 0;
              0, cos(theta2), sin(theta2);
              0, -sin(theta2), cos(theta2)];
    elseif(axis(2) == 2)
        R2 = [cos(theta2), 0, -sin(theta2);
              0, 1, 0;
              sin(theta2), 0, cos(theta2)];
    elseif(axis(2) == 3)
        R2 = [cos(theta2), sin(theta2), 0;                       
              -sin(theta2), cos(theta2), 0;
              0, 0, 1];
    end

    % theta 3.
    if(axis(3) == 1)
        R3 = [1, 0, 0;
              0, cos(theta3), sin(theta3);
              0, -sin(theta3), cos(theta3)];
    elseif(axis(3) == 2)
        R3 = [cos(theta3), 0, -sin(theta3);
              0, 1, 0;
              sin(theta3), 0, cos(theta3)];
    elseif(axis(3) == 3)
        R3 = [cos(theta3), sin(theta3), 0;                       
              -sin(theta3), cos(theta3), 0;
              0, 0, 1];
    end

    % compute R
    R = R3 * R2 * R1;

end