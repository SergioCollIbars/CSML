function [STM] = integrateSTM(TIME, state, asterParams, poleParams, Cnm, Snm, order)
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
%           r_ACI: position vector in the inertial frame                  %
%           planetParams: planet parameters                               %
%                     [GM, Re, nmax, normalized]                          %
%           poleParams: pole parameters                                   %
%                   [W, W0, RA, DEC]                                      %
%           Cmat: SH C coefficients                                       %
%           Smat: SM S coefficients                                       % 
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
        % construct Jacobian
        [ddU_ACI, ~, ~] = gradiometer_meas(TIME(j) ,asterParams, poleParams, state(j, :), ...
                zeros(9, 1), Cnm, Snm, eye(3,3));
        T = [ddU_ACI(1), ddU_ACI(2), ddU_ACI(3);...
            ddU_ACI(4), ddU_ACI(5), ddU_ACI(6);...
            ddU_ACI(7), ddU_ACI(8), ddU_ACI(9)];
        
        J = [zeros(3, 3), eye(3,3); T, zeros(3, 3)];

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

