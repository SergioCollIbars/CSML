function [Y] = add_angularComponents(Y, attitude, At, angVel, angAcc)
    % Description: Given the gradiometer measurements in the ACI frame,
    % rotate them to the body frame and include the angular velocity and
    % acceleration Dyads to the measurements.
    Nt  = length(Y(1, :));
    for j = 1:Nt
      omega = [0, angVel(3, j), -angVel(2, j);...
        -angVel(3, j), 0 ,angVel(1, j);...
         angVel(2, j), -angVel(1, j), 0];
      omegaDot = [0, angAcc(3, j), -angAcc(2, j);...
        -angAcc(3, j), 0 ,angAcc(1, j);...
        angAcc(2, j), -angAcc(1, j), 0];
      a = omega^2;
      b = omegaDot;
      A = [a(1,1); a(1,2); a(1,3); a(2, 1); a(2, 2); a(2, 3); ...
          a(3, 1); a(3, 2); a(3, 3)];
      B = [b(1,1); b(1,2); b(1,3); b(2, 1); b(2, 2); b(2, 3); ...
          b(3, 1); b(3, 2); b(3, 3)];

      % rotate to actual body frame. 
      yaw  = attitude(1, j); pitch = attitude(2,j); roll = attitude(3, j);
      BN =rotationMatrix(yaw, pitch, roll, ...
          [3, 2, 1]);   % from ACI to Nominal body frame
      AB =rotationMatrix(At(1, j), At(2, j), At(3, j), ...
          [3, 2, 1]);   % from Nominal body frame to actual frame (A)
      AN = AB * BN;

      jN = [Y(1, j), Y(2, j), Y(3, j);Y(4, j), Y(5, j), Y(6, j);...
          Y(7, j), Y(8, j), Y(9, j)];
      jB = AN * jN * AN';
      JB = [jB(1,1); jB(1,2) ; jB(1,3); jB(2,1);...
         jB(2,2); jB(2,3) ; jB(3,1); jB(3,2); jB(3,3)];
      Y(:, j) = JB + A - B;
    end
end

