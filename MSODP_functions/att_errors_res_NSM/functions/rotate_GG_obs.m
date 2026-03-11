function [GG_rot, GG_org_cmpl] = rotate_GG_obs(GG_org,orientation)
    % Rotate gradiometer observations given the set of 3-2-1 rotation
    % angles.
    % input GG_org (Nt x 7): time, xx, xy, xz, yy, yz, zz
    %       orientation (Nt x 3): yaw, pitch, roll
    % output GG_rot (Nt x 10): time, xx, xy, xz, yx, yy, yz, zx, zy, zz
    %        GG_org_cmpl (Nt x 10): time, xx, xy, xz, yx, yy, yz, zx, zy, zz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Nt = length(GG_org(:, 1));
    GG_rot      = nan(Nt, 10); GG_rot(:, 1) = GG_org(:, 1);
    GG_org_cmpl = GG_rot;
    for k = 1:Nt
        yaw = orientation(1, k); pitch = orientation(2, k);
        roll= orientation(3, k);

        TN   = [GG_org(k, 2), GG_org(k, 3), GG_org(k, 4);...
             GG_org(k, 3),  GG_org(k, 5),  GG_org(k, 6);...
             GG_org(k, 4),  GG_org(k, 6),  GG_org(k, 7)];

        BN = rotationMatrix(yaw, pitch, roll, [3, 2, 1]);
        TB = BN  * TN * BN';

        GG_rot(k, 2:end) = ...
            [TB(1,1);TB(1,2);TB(1,3);TB(2,1);TB(2,2);TB(2,3);...
             TB(3,1);TB(3,2);TB(3,3)];

        GG_org_cmpl(k, 2:end) = ...
            [TN(1,1);TN(1,2);TN(1,3);TN(2,1);TN(2,2);TN(2,3);...
             TN(3,1);TN(3,2);TN(3,3)];
    end

end

