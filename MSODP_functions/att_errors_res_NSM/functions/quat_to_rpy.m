function [roll, pitch, yaw] = quat_to_rpy(q)
%QUAT_TO_RPY_SCALAR_FIRST Convert quaternion [q0 q1 q2 q3] to roll/pitch/yaw
%
% Convention:
%   q = [q0 qx qy qz], scalar first
%   Returns 3-2-1 Euler angles:
%       roll  = rotation about x
%       pitch = rotation about y
%       yaw   = rotation about z
%
% Angles in radians.

    q0 = q(:,1);
    qx = q(:,2);
    qy = q(:,3);
    qz = q(:,4);

    % Roll
    sinr_cosp = 2 .* (q0 .* qx + qy .* qz);
    cosr_cosp = 1 - 2 .* (qx.^2 + qy.^2);
    roll = wrapTo2Pi(atan2(sinr_cosp, cosr_cosp));

    % Pitch
    sinp = 2 .* (q0 .* qy - qz .* qx);
    sinp = max(min(sinp, 1), -1);   % protect against roundoff
    pitch = asin(sinp);

    % Yaw
    siny_cosp = 2 .* (q0 .* qz + qx .* qy);
    cosy_cosp = 1 - 2 .* (qy.^2 + qz.^2);
    yaw = wrapTo2Pi(atan2(siny_cosp, cosy_cosp));
end