function [r_ACI, v_ACI] = Kepler_motion(alpha0, mu, t)
    % ------------------------------------------------------------------- %
    %                         PACKAGE FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 19/09/2022

    % Description: package to compute the orbit at each time.

    % Input:
    %   t: current time vector
    %   mu: gravitational parameter
    %   tau: periapsis pass time
    %   alpha0: initial orbital elemets

    % Output:
    %   r_ACI: new position vector in {I, J, K} frame
    %   v_ACI: new velocity vector in {I, J, K} frame

    % ------------------------------------------------------------------- %
    
    % Extrac inital orbital elements
    Omega = alpha0(4);
    i = alpha0(3);
    omega = alpha0(5);
    e = alpha0(2);
    a = alpha0(1);
    rho = a * (1 - e^2);
    f0 = alpha0(6);

    % Rotational matrix {e, ep} -> {I, J, K}
    R_Omega = [cos(Omega) sin(Omega) 0;...
               -sin(Omega) cos(Omega) 0;...
               0 0 1];
    R_i = [1 0 0;...
           0 cos(i) sin(i);...
           0 -sin(i) cos(i)];
    R_omega = [cos(omega) sin(omega) 0;...
               -sin(omega) cos(omega) 0;...
               0 0 1];
    BN = R_omega * R_i * R_Omega;

    % compute mean motion
    n = sqrt(mu / a^3);

    % time through periapsis
    sE0 = sqrt(1 - e^2) * sin(f0) / (1 + e*cos(f0));
    cE0 = (e + cos(f0)) / (1 + e*cos(f0));
    E0 = atan2(sE0, cE0);
    
    tau = - 1/n * (E0 - e*sin(E0));

    % Compute Mean anomaly.
    M = n * (t - tau);

    % Solve Kepler equation
    E = keplerEq_solver(M, e);

    % Compute true anomaly.
    s_f = sqrt(1 - norm(e)^2) * sin(E) / (1 - norm(e) * cos(E));
    c_f = (cos(E) - norm(e)) / (1 - norm(e) * cos(E));

    f = atan2(s_f, c_f);                                % New true anomaly.

    % compute new values of position and velocity. frame = {e, ep}
    r = rho / (1 + norm(e) * cos(f)) * [cos(f); sin(f); 0];          
    v = sqrt(mu /rho) * [- sin(f); norm(e) + cos(f); 0];    

    % convert modules to {I, J, K} frame of reference.
    r_ACI = BN' * r;
    v_ACI = BN' * v;
end

%%                          FUNCTIONS

function [E] = keplerEq_solver(M, e)
    % ------------------------------------------------------------------- %
    %                          KEPLER EQ FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 19/09/2022

    % Description: Function to compute the Kepler equation using Newtons 
    %   method.

    % Input:
    %   M: Mean anomaly
    %   e: eccentricity vector

    % Output:
    %   E: Eccentric anomaly
    % ------------------------------------------------------------------- %
    
    Delta = 1E-12;                       % Convergence value.
    error = Delta + 1;                  % Function error. Initial value.

    k = 1;                              % Iteration value.
    while(abs(error) > Delta)
        
        if (k == 1)
            E = M;
        else
            E = E + error;
        end

        F =  M - (E  - norm(e) * sin(E));      % function value.

        dE = - 1 + norm(e) * cos(E);
        error = - F / dE;

        k = k + 1;                      % Update iteration value.
    end
end