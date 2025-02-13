function [rFix,vFix] = rot2iner(th,rRot,vRot,mu,n,a,ref)
% Transformation from M body-centered rotating (nondim.) coordinates to inertial (nondim.) coordinates
% Input:    th      - vector    (rotation angle t*n)
% Input:    rRot    - matrix    (3xN matrix with position in rotating frame)
% Input:    vRot    - matrix    (3xN matrix with velocity in rotating frame)
% Input:    ref     - scalar    (objective reference frame: 0 for barycentric, 1 for primary mass centered)
% Oitput:   rFix    - matrix    (3xN matrix with position in inertial frame)
% Output:   vFix    - matrix    (3xN matrix with velocity in inertial frame)
% 
% Elias Marin (revised 23-02-2024)

%% Initial settings
if ref == 0
    dist = 0;
elseif ref == 1
    dist = -mu* a;
end

%% Computations
rFix = zeros(size(rRot));
vFix = zeros(size(vRot));

% Rotation angles
c = cos(th);
s = sin(th);

% Position rotation
rFix(1,:) =  c.*(rRot(1,:) - dist*n) - s.*rRot(2,:);
rFix(2,:) = s.*(rRot(1,:) - dist*n) + c.*rRot(2,:);
rFix(3,:) = rRot(3,:);

% Velocity rotation
vFix(1,:) =  c.*(vRot(1,:) - rRot(2,:)*n) + s.*(vRot(2,:)+rRot(1,:)*n - dist*n);
vFix(2,:) = -s.*(vRot(1,:) - rRot(2,:)*n) + c.*(vRot(2,:)+rRot(1,:)*n - dist*n);
vFix(3,:) = vRot(3,:);


end
