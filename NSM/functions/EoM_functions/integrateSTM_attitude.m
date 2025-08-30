function [STM] = integrateSTM_attitude(TIME, state, order, Iner)
%%                      RECONSTRUCT TRAJECTORY                           %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/31/2024                                                      %
%                                                                         %
%   Description: function to integrate gradiometer measurements and       %
%       reconstruct relative trajectory.                                  %
%                                                                         %
%   Inputs:                                                               %
%           TIME: time vector to evaluate meas                            %
%           state: Euler anlges + angular velocity                        %
%           Iner: Inertia matrix                                          %
%                                                                         %     
%   Output: STM:  in ACI frame                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nt = length(TIME);

    % STM definition
    PHI = eye(6, 6);

    % computed STM
    STM = ones(Nt, 36) * NaN;
    STM(1, :) = reshape(PHI, [1, 36]);

    % Adam bashforth initialization
    x = 0;
    f = 0;
    
    % loop
    for j = 1:Nt-1
        theta = state(2, j);
        phi   = state(3, j);
        omega = state(4:6, j);
        w1 = omega(1); w2 = omega(2); w3 = omega(3);

        A  = [-sin(theta), 0, 1;...
                sin(phi)*cos(theta), cos(phi), 0;...
                cos(phi)*cos(theta), -sin(phi), 0];
        qDot_wrt_omega = inv(A);
    
        [t1] = compute_skewMat(omega);
        [t2] = compute_skewMat(Iner * omega);
        omegaDot_wrt_omega = -inv(Iner) * (t1 * Iner - t2);

        PsiDot   = 1/cos(theta) * (sin(phi)*w2 + cos(phi)*w3);
        ThetaDot = 1/cos(theta) * (cos(phi)*cos(theta)*w2 - sin(phi)*cos(theta)*w3);
        PhiDot   = 1/cos(theta) * (w1 + sin(phi)*sin(theta)*w2 + cos(phi)*sin(theta)*w3);
    
        a = sec(theta)*tan(theta)*ThetaDot*cos(theta) + ...
            1/cos(theta)*(-cos(phi)*sin(theta)*w2 + sin(phi)*sin(theta)*w3);
        b = sec(theta)*tan(theta)*PhiDot*cos(theta) + ...
            1/cos(theta)*(sin(phi)*cos(theta)*w2 + cos(phi)*cos(theta)*w3);
    
        C = [0, sec(theta)*tan(theta)*PsiDot*cos(theta), 1/cos(theta)*(cos(phi)*w2 - sin(phi)*w3);...
             0, a, 1/cos(theta)*(-sin(phi)*cos(theta)*w2 - cos(phi)*cos(theta)*w3);...
             0, b, 1/cos(theta)*(cos(phi)*sin(theta)*w2 - sin(phi)*sin(theta)*w3)];

        % construct Jacobian
        J = [C, qDot_wrt_omega; zeros(3, 3), omegaDot_wrt_omega];

        % integrate PHI
        t_span = [TIME(j), TIME(j+1)];
        At = diff(t_span);

        PHI_dot = J * PHI;

        % integrate STM
        if(order == 0)
            [PHI_new, PHI] = integration_order0(PHI, PHI_dot, At);
        elseif(order == 2)
            [PHI_new, PHI, x, f] = integration_order2(x, f, PHI, PHI_dot, ...
                At, j);
        elseif(order == 3)
            [PHI_new, PHI, x, f] = integration_order3(x, f, PHI, PHI_dot, ...
                At, j);
        elseif(order == 4)
            [PHI_new, PHI, x, f] = integration_order4(x, f, PHI, PHI_dot, ...
                At, j);
        elseif(order == 5)
            [PHI_new, PHI, x, f] = integration_order5(x, f, PHI, PHI_dot, ...
                At, j);
        end
        
        % store STM
        STM(j+1, :) = reshape(PHI_new, [1, 36]);
    end
end

%% functions

function [PHI_new, PHI] = integration_order0(PHI, PHI_dot, At)
    [PHI_new] = euler_integration(PHI, PHI_dot, At);
    PHI = PHI_new;
end

function [PHI_new, PHI, x, f] = integration_order2(x, f, PHI, PHI_dot, ...
    At, j)
        if(j == 1)
            % Euler integration method
            [PHI_new] = euler_integration(PHI, PHI_dot, At);
        else
            % Adams Bashforth integration method
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 2);
        end
        x = reshape(PHI_new, [36,1]);
        f = reshape(PHI_dot, [36, 1]);
        PHI = PHI_new;
end

function [PHI_new, PHI, x, f] = integration_order3(x, f, PHI, PHI_dot, ...
    At, j)
        if(j == 1)
            % Euler integration method
            [PHI_new] = euler_integration(PHI, PHI_dot, At);
            f = reshape(PHI_dot, [36, 1]);
        elseif(j == 2)
            % Adams Bashforth integration method. Order 2
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 2);
            f = [reshape(PHI_dot, [36, 1]); f];
        else
            % Adams Bashforth integration method. Order 3
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 3);
            f = [reshape(PHI_dot, [36, 1]); f(1:36)];
        end
        x = reshape(PHI_new, [36,1]);
        PHI = PHI_new;
end

function [PHI_new, PHI, x, f] = integration_order4(x, f, PHI, PHI_dot, ...
    At, j)
        if(j == 1)
            % Euler integration method
            [PHI_new] = euler_integration(PHI, PHI_dot, At);
            f = reshape(PHI_dot, [36, 1]);
        elseif(j == 2)
            % Adams Bashforth integration method. Order 2
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 2);
            f = [reshape(PHI_dot, [36, 1]); f];
        elseif(j == 3)
            % Adams Bashforth integration method. Order 3
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 3);
            f = [reshape(PHI_dot, [36, 1]); f];
        else
            % Adams Bashforth integration method. Order 4
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 4);
            f = [reshape(PHI_dot, [36, 1]); f(1:72)];
        end
        x = reshape(PHI_new, [36,1]);
        PHI = PHI_new;
end

function [PHI_new, PHI, x, f] = integration_order5(x, f, PHI, PHI_dot, ...
    At, j)
        if(j == 1)
            % Euler integration method
            [PHI_new] = euler_integration(PHI, PHI_dot, At);
            f = reshape(PHI_dot, [36, 1]);
        elseif(j == 2)
            % Adams Bashforth integration method. Order 2
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 2);
            f = [reshape(PHI_dot, [36, 1]); f];
        elseif(j == 3)
            % Adams Bashforth integration method. Order 3
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 3);
            f = [reshape(PHI_dot, [36, 1]); f];
        elseif(j == 4)
            % Adams Bashforth integration method. Order 4
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 4);
            f = [reshape(PHI_dot, [36, 1]); f];
        else
            % Adams Bashforth integration method. Order 5
            [PHI_new] = adamsBashforth_integration(x, [reshape(PHI_dot, [36, 1]); f], ...
                At, 5);
            f = [reshape(PHI_dot, [36, 1]); f(1:108)];
        end
        x = reshape(PHI_new, [36,1]);
        PHI = PHI_new;
end


function [skw] = compute_skewMat(vec)
    skw = [0, -vec(3), vec(2);...
           vec(3), 0, -vec(1);...
          -vec(2), vec(1), 0];
end
