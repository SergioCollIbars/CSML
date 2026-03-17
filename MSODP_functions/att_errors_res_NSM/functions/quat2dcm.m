function C = quat2dcm(q)
%QUAT2DCM_CUSTOM Convert quaternion to DCM
%
% Input:
%   q : 4x1 or 1x4 quaternion [q0 q1 q2 q3]
%       where q0 is the scalar part
%
% Output:
%   C : 3x3 direction cosine matrix
%
% Note:
%   The quaternion is normalized internally.

    q = q(:);  % force column vector

    if length(q) ~= 4
        error('Input quaternion must have 4 elements.');
    end

    % Normalize quaternion
    q = q / norm(q);

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    C = [ 1 - 2*(q2^2 + q3^2),   2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2);
          2*(q1*q2 + q0*q3),     1 - 2*(q1^2 + q3^2),   2*(q2*q3 - q0*q1);
          2*(q1*q3 - q0*q2),     2*(q2*q3 + q0*q1),     1 - 2*(q1^2 + q2^2) ];
end