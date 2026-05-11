function [H_omegaDot] = angularAcc_dyad_partials(dt, omega, I)
    
    d_dwx = (1/dt)*[0,0,0,0,0,-1,0,1,0]'; 
    d_dwy = (1/dt)*[0,0,1,0,0,0,-1,0,0]'; 
    d_dwz = (1/dt)*[0,-1,0,1,0,0,0,0,0]'; 
    
    H_omega_stack = [d_dwx, d_dwy, d_dwz];


    Iner = I;    % inertia matrix
    
    [t1] = compute_skewMat(omega);
    [t2] = compute_skewMat(Iner * omega);
    H_omegaDot = -inv(Iner) * (t1 * Iner - t2);

    H_omegaDot = [zeros(1, 3); - H_omegaDot(3, :); H_omegaDot(2, :); ...
        H_omegaDot(3, :); zeros(1, 3); - H_omegaDot(1, :); ...
        -H_omegaDot(2, :); H_omegaDot(1, :); zeros(1, 3)]; % Dyad partials w.r.t angular velocity
end

function [skw] = compute_skewMat(vec)
    skw = [0, -vec(3), vec(2);...
           vec(3), 0, -vec(1);...
          -vec(2), vec(1), 0];
end