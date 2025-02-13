function [E, E0, phiE, phiE_left, phiE_up] = ...
    load_extraParams(mode, Nc, Ns)
    
    % mode 0. No extra parameters
    if(mode == 0)
        E = zeros(10, 1);
        E0 = zeros(10, 1);
        phiE = 0;
        phiE_left = 0;
        phiE_up = 0;
    end
    % mode 1. Noise parameters only (CKF)
    if(mode == 1)
        E = zeros(2, 1);
        d = 1E-16;
        % phase A values
        E0 = [-1*1E-9/2, d, 1*1E-9/2, d, 1*1E-9/2, 1E-3 * d, ...
            -1*1E-9/2, d, -1*1E-9/2, d, -1*1E-9/2, d]';
        
% % % %         % phase B values
% %         E0 = [-1.41E-7, 3.88E-17, 4.47E-9, 5.23E-18, 8.21E-9, 1.32E-19, ...
% %             -3.17E-7, 1.12E-15, -9.75E-6, 5.27E-15, -5.65E-8, 1.88E-17]';

% %         % super close values
% %         E0 = [-282*1E-9/2, 4E-17, 8.88*1E-9/2, 6E-17, 16.42*1E-9/2, 6E-19, ...
% %             -636*1E-9/2, 1E-15, -19500*1E-9/2, 6E-15,-113*1E-9/2, 2E-17]'; 

        % Jacobian and STM dot matrices
        Ae = [0, 1;0, 0];
        phiE = blkdiag(Ae, Ae, Ae, Ae, Ae, Ae);
        phiE_left = zeros(12, Nc + Ns);
        phiE_up = phiE_left';
    end
    % mode 2. Noise parameters withou apriori (LS)
    if(mode == 2)
        E  = zeros(12, 1);
        E0 = ones(12, 1);
        phiE = 0;
        phiE_left = 0;
        phiE_up = 0;
    end
    
end

