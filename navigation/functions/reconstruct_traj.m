function [X_new, STM] = reconstruct_traj(TIME, ddU_ACI, X0,...
    enviroment, order)
%%                      RECONSTRUCT TRAJECTORY                           %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/31/2024                                                      %
%                                                                         %
%   Description: function to integrate gradiometer measurements and       %
%       reconstruct relative trajectory.                                  %
%                                                                         %
%   Inputs: TIME: time vector                                             %
%           r_ACI: position vector in the inertial frame                  %
%           planetParams: planet parameters                               %
%                     [GM, Re, nmax, normalized]                          %
%           poleParams: pole parameters                                   %
%                   [W, W0, RA, DEC]                                      %
%           Cmat: SH C coefficients                                       %
%           Smat: SM S coefficients                                       % 
%                                                                         %     
%   Output: ddU_ACI: gradiometer measurement in ACI frame                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nt = length(TIME);

    % STM definition
    PHI = eye(6, 6);

    % reconstructed trajectory 
    X_new = ones(6, Nt) * NaN;
    X_new(:, 1) = X0;
    
    % computed STM
    STM = ones(Nt, 36) * NaN;
    STM(1, :) = reshape(PHI, [1, 36]);

    % Adam bashforth initialization
    x = 0;
    f = 0;
    
    % loop
    for j = 1:Nt-1
        % construct Jacobian
        T = reshape(ddU_ACI(:, j), [3,3]);
        if(enviroment == "CR3BP")
            J = [zeros(3, 3), eye(3,3); T, [0,2,0;-2,0,0;0,0,0]];
        elseif(enviroment == "2BP" || enviroment == "CR3BP_inertial")
            J = [zeros(3, 3), eye(3,3); T, zeros(3, 3)];
        end

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

        % update state
        X_new(:, j+1) = PHI_new  * X_new(:, 1);

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

