function [phi, phiDot, phiDdot, theta,thetaDot, thetaDdot, psi, psiDot, psiDdot] = attitude_profile(PM, axis, t, val)
    % ------------------------------------------------------------------- %
    %                     ATTIUDE PROFILE FUNCTION
    %
    %   Input:
    %       PM: profile mode parameter
    %       axis: array with the axis affected
    %       t: time vector
    % ------------------------------------------------------------------- %
    
    % get axis affected
    axisN = length(axis);
    axis1 = 0;
    axis2 = 0;
    axis3 = 0;

    for i = 1: axisN
        if(axis(i) == 1)
            axis1 = 1;
        elseif(axis(i) == 2)
            axis2 = 1;
        elseif(axis(i) == 3)
            axis3 = 1;
        end
    end

    % Euler angles matrix definition
    phi = zeros(1, length(t));          % R[AXIS 1]
    theta = zeros(1, length(t));        % R[AXIS 2]
    psi = zeros(1, length(t));          % R[AXIS 3]
    
    phiDot = zeros(1, length(t));       % R[AXIS 1]
    thetaDot = zeros(1, length(t));     % R[AXIS 2]
    psiDot = zeros(1, length(t));       % R[AXIS 3]
    
    phiDdot = zeros(1, length(t));      % R[AXIS 1]
    thetaDdot = zeros(1, length(t));    % R[AXIS 2]
    psiDdot = zeros(1, length(t));      % R[AXIS 3]
    
    % compute attitude profile
    if (PM == "constant")
        
        for k =1:length(t)
            if(axis1 == 1)
                phi(k)      = val;
                phiDot(k)   = 0;
                phiDdot(k)  = 0;
            end
            if(axis2 == 1)
                theta(k)    = val;
                thetaDot(k) = 0;
                thetaDdot(k)= 0;
            end
            if(axis3 == 1)
                psi(k)      = val;
                psiDot(k)   = 0;
                psiDdot(k)  = 0;
            end
        end
        
     elseif(PM == "linear")
        shift = 0;

        for k =1:length(t)
            if(axis1 == 1)
                phi(k)      = val * t(k) + shift;
                phiDot(k)   = val;
                phiDdot(k)  = 0;
            end
            if(axis2 == 1)
                theta(k)    = val * t(k) + shift;
                thetaDot(k) = val;
                thetaDdot(k)= 0;
            end
            if(axis3 == 1)
                psi(k)      = val * t(k) + shift;
                psiDot(k)   = val;
                psiDdot(k)  = 0;
            end
        end

    elseif(PM == "null")
        
        for k =1:length(t)
            if(axis1 == 1)
                phi(k)      = 0;
                phiDot(k)   = 0;
                phiDdot(k)  = 0;
            end
            if(axis2 == 1)
                theta(k)    = 0;
                thetaDot(k) = 0;
                thetaDdot(k)= 0;
            end
            if(axis3 == 1)
                psi(k)      = 0;
                psiDot(k)   = 0;
                psiDdot(k)  = 0;
            end
        end
    elseif(PM == "sinusoidal")
        
         for k =1:length(t)
            if(axis1 == 1)
                phi(k)      = val * sin(t(k)); % add orbit rate
                phiDot(k)   = val * cos(t(k));
                phiDdot(k)  = -val * sin(t(k));
            end
            if(axis2 == 1)
                theta(k)    = val * sin(t(k));
                thetaDot(k) = val * cos(t(k));
                thetaDdot(k)= -val * sin(t(k));
            end
            if(axis3 == 1)
                psi(k)      = val * sin(t(k));
                psiDot(k)   = val * cos(t(k));
                psiDdot(k)  = -val * sin(t(k));
            end
        end
    elseif(PM == "WGN")     % HOW TO MODEL THE DERIVATIVE? Numerically?
        for k =1:length(t)
            if(axis1 == 1)
                phi(k)      = wgn(1,1,0);
                phiDot(k)   = val * cos(t(k));
                phiDdot(k)  = -val * sin(t(k));
            end
            if(axis2 == 1)
                theta(k)    = val * sin(t(k));
                thetaDot(k) = val * cos(t(k));
                thetaDdot(k)= -val * sin(t(k));
            end
            if(axis3 == 1)
                psi(k)      = val * sin(t(k));
                psiDot(k)   = val * cos(t(k));
                psiDdot(k)  = -val * sin(t(k));
            end
        end
    end

end