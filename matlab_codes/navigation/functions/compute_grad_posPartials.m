function [ddU_dxyz] = compute_grad_posPartials(GM, x, y, z)
    %%      COMPUTE GRADIOMETER POSITION PARTIALS
    %   Description: compute gradiometer position partials in the inertial
    %   frame.
    %   Input: GM = point mass parameter
    %          x, y, z = position coordinates in the inertial frame
    %   Output: ddU_dxyz = grad partials w.r.t position.
    %           ddU_dxyz = [dU11_dx, dU11_dy, dU11_dz; ...
    %                       dU12_dx, dU12_dy, dU12_dz; ... etc]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % position radius norm
    r = sqrt(x*x + y*y + z*z);

    % radius power to 5 and power to 7
    r5 = r^5; r7 = r^7;

    % gradiometer partials
    dU11_dx = -9*GM/r5*x + 15*GM/r7 *x^3;
    dU11_dy = -3*GM/r5*y + 15*GM/r7*x^2*y;
    dU11_dz = -3*GM/r5*z + 15*GM/r7*x^2*z;

    dU12_dx = -3*GM/r5*y + 15*GM/r7 *x^2*y;
    dU12_dy = -3*GM/r5*x + 15*GM/r7*y^2*x;
    dU12_dz = 15*GM/r7*x*y*z;

    dU13_dx = -3*GM/r5*z + 15*GM/r7 *x^2*z;
    dU13_dy = 15*GM/r7*x*y*z;
    dU13_dz = -3*GM/r5*x + 15*GM/r7*z^2*x;

    dU22_dy = -9*GM/r5*y + 15*GM/r7 *y^3;
    dU22_dx = -3*GM/r5*x + 15*GM/r7*y^2*x;
    dU22_dz = -3*GM/r5*z + 15*GM/r7*y^2*z;

    dU23_dx = 15*GM/r7*x*y*z;
    dU23_dy = -3*GM/r5*z + 15*GM/r7*y^2*z;
    dU23_dz = -3*GM/r5*y + 15*GM/r7*z^2*y;


    dU33_dx = -3*GM/r5*x + 15*GM/r7*z^2*x;
    dU33_dy = -3*GM/r5*y + 15*GM/r7*z^2*y;
    dU33_dz = -9*GM/r5*z + 15*GM/r7 *z^3;

    % output
    ddU_dxyz = -[dU11_dx, dU11_dy, dU11_dz;...
        dU12_dx, dU12_dy, dU12_dz;...
        dU13_dx, dU13_dy, dU13_dz;...
        dU22_dx, dU22_dy, dU22_dz;...
        dU23_dx, dU23_dy, dU23_dz;...
        dU33_dx, dU33_dy, dU33_dz];
end

