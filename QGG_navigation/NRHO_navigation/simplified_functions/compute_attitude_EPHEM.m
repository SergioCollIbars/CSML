function [BN_matrix] = compute_attitude_EPHEM(t, frec, Nn, attitude)
    % Compute attitude using the EPHEMERIDES model.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % time points
    Nt = length(t);
  
    % compute Euler angles noise
    sigmaAng = [1;1;1] * 1 * sqrt(frec) * pi / (180 * 3600);        % [rad]

    noise_ST = ones(3, Nt) * NaN;
    for k = 1:length(sigmaAng)
         noise_ST(k, :) = normrnd(0, sigmaAng(k), [1, Nt]) * Nn;
    end

    % compute attitude inertial (yaw - pitch - roll)
    Angles = zeros(3, Nt);
    if(attitude == "RTN")
        Angles(1, :) = ones(1, Nt).*pi/4;
        Angles(2, :) = ones(1, Nt).*pi/6;
        Angles(3, :) = ones(1, Nt).*pi/4;
    end
    
    BN_matrix = ones(3 * Nt, 3) * NaN;
    for k = 1:Nt
        Psi   = Angles(1, k) + noise_ST(1, k) * Nn;   % yaw
        theta = Angles(2, k) + noise_ST(2, k) * Nn;   % pitch
        phi   = Angles(3, k) + noise_ST(3, k) * Nn;   % roll

        BN = [cos(theta)*cos(Psi),-cos(theta)*sin(Psi),sin(theta);...
            sin(phi)*sin(theta)*cos(Psi)+cos(phi)*sin(Psi), ...
            -sin(Psi)*sin(theta)*sin(Psi)+cos(phi)*cos(Psi), -sin(phi)*cos(theta);...
            sin(phi)*sin(Psi)-cos(phi)*cos(Psi)*sin(theta), ...
            sin(phi)*cos(Psi)+cos(phi)*sin(theta)*sin(Psi),...
            cos(phi)*cos(theta)];
        
        maxInd = 3 * k; minInd = maxInd - 2;
        BN_matrix(minInd:maxInd, :) = BN;
    end
end

