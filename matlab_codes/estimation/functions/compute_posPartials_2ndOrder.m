function [ddU_dxyz] = compute_posPartials_2ndOrder(GM, x, y, z)
    %%      COMPUTE GRADIOMETER POSITION PARTIALS
    %   Description: compute gradiometer 2nd order position partials 
    %           in the inertial frame.
    %   Input: GM = point mass parameter
    %          x, y, z = position coordinates in the inertial frame
    %   Output: ddU_dxyz = grad partials w.r.t position.
    %           ddU_dxyz = [ddU11_dxx, ddU11_dxy, ddU11_dxz, ddU11_dyy, ddU11_dyz, ddU11_dzz; ...
    %                       ddU12_dxx, ddU12_dxy, ddU12_dxz, ddU12_dyy, ddU12_dyz, ddU12_dzz; ... etc]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % vecnor norm
    r = sqrt(x^2 + y^2 + z^2);
    r5 = r^5;
    r7 = r^7;
    r9 = r^9;

    % U11 2nd order partials
    ddU11_dxx = -9*GM/r5 + 90*GM/r7*x^2 - 105*GM/r9*x^4;
    ddU11_dyy = -3*GM/r5 + 15*GM/r7*(x^2 + y^2) - 105*GM/r9*x^2*y^2;
    ddU11_dzz = -3*GM/r5 + 15*GM/r7*(x^2 + z^2) -105*GM/r9*x^2*z^2;
    ddU11_dxy = 45*GM/r7*x*y - 105*GM/r9*x^3*y;
    ddU11_dxz = 45*GM/r7*x*z - 105*GM/r9*x^3*z;
    ddU11_dyz = 15*GM/r7*y*z -105*GM/r9*x^2*y*z;

    % U12 2nd order partials
    ddU12_dxx = 45*GM/r7*x*y -105*GM/r9*x^3*y;
    ddU12_dyy = 45*GM/r7*x*y - 105*GM/r9*y^3*x;
    ddU12_dzz = 15*GM/r7*x*y - 105*GM/r9*x*y*z^2;
    ddU12_dxy = -3*GM/r5 + 15*GM/r7*(y^2 + x^2) -105*GM/r9*x^2*y^2;
    ddU12_dxz = 15*GM/r7*y*z - 105*GM/r9*x^2*y*z;
    ddU12_dyz = 15*GM/r7*x*z - 105*GM/r9*y^2*x*z;

    % U13 2nd order partials
    ddU13_dxx = 45*GM/r7*x*z -105*GM/r9*x^3*z;
    ddU13_dyy = 15*GM/r7*x*z -105*GM/r9*x*z*y^2;
    ddU13_dzz = 45*GM/r7*x*z -105*GM/r9*x^2*z^2;
    ddU13_dxy = 15*GM/r7*z*y - 105*GM/r9*x^2*y*z;
    ddU13_dxz = -3*GM/r5 + 15*GM/r7*(z^2+x^2) -105*GM/r9*x^2*z^2;
    ddU13_dyz = 15*GM/r7*x*y -105*GM/r9*x*y*z^2;

    % U22 2nd order partials
    ddU22_dxx = -3*GM/r5 + 15*GM/r7*(x^2 + y^2) -105*GM/r9*x^2*y^2;
    ddU22_dyy = -9*GM/r5 + 90*GM/r7*y^2 -105*GM/r9*y^4;
    ddU22_dzz = -3*GM/r5 + 15*GM/r7*(z^2+y^2) -105*GM/r9*y^2*z^2;
    ddU22_dxy = 45*GM/r7*x*y -105*GM/r9*y^3*x;
    ddU22_dxz = 15*GM/r7*x*z -105*GM/r9*y^2*x*z;
    ddU22_dyz = 45*GM/r7*y*z -105*GM/r9*y^3*z;

    % U23 2nd order partials
    ddU23_dxx = 15*GM/r7*y*z -105*GM/r9*x^2*y*z;
    ddU23_dyy = 45*GM/r7*z*y -105*GM/r9*y^3*z;
    ddU23_dzz = 45*GM/r7*z*y -105*GM/r9*z^3*y;
    ddU23_dxy = 15*GM/r7*x*z -105*GM/r9*x*y^2*z;
    ddU23_dxz = 15*GM/r7*x*y -105*GM/r9*x*y*z^2;
    ddU23_dyz = -3*GM/r5 + 15*GM/r7*(z^2 + y^2) -105*GM/r9*y^2*z^2;

    % ensamble matrix
    ddU_dxyz = [0.5*ddU11_dxx, ddU11_dxy, ddU11_dxz, 0.5*ddU11_dyy, ddU11_dyz, 0.5*ddU11_dzz;...
        0.5*ddU12_dxx, ddU12_dxy, ddU12_dxz, 0.5*ddU12_dyy, ddU12_dyz, 0.5*ddU12_dzz;...
        0.5*ddU13_dxx, ddU13_dxy, ddU13_dxz, 0.5*ddU13_dyy, ddU13_dyz, 0.5*ddU13_dzz;...
        0.5*ddU22_dxx, ddU22_dxy, ddU22_dxz, 0.5*ddU22_dyy, ddU22_dyz, 0.5*ddU22_dzz;...
        0.5*ddU23_dxx, ddU23_dxy, ddU23_dxz, 0.5*ddU23_dyy, ddU23_dyz, 0.5*ddU23_dzz];
end

