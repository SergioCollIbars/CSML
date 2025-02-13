function [dx] = dynamic_propagator(t, X, varObj)
    %%                   DYNAMICS PROPAGATOR FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 01/20/2023
    %
    %   Description: This function gives the diferential equation for the
    %   state and STM propagation and dynamics.
    %
    %   Input:
    %       t: time vector
    %       X: state
    %  
    %   Output: 
    %       dx: diferential equation set
    % --------------------------------------------------------------------%
    
    % extract values
    Np = varObj.Np;
    Nc = varObj.Nc;
    Ns = varObj.Ns;

    % extract STM
    PHI = reshape(X(Np+1:end), [Np, Np]);
    
    % ensamble Jacobian matrix
    if(Np == (Nc + Ns))     % no extra parameters
        A = zeros(Np, Np);
    else    % there's extra parameters
        % get extra parameters
        [~, ~, phiE, phiE_left, phiE_up] = ...
            load_extraParams(varObj.extraParam, Nc, Ns);
       
        A = [zeros(Nc+Ns, Nc+Ns), phiE_up;...
            phiE_left, phiE];
    end
    N_phi = length(A(:, 1));

    % update PHI
    PHI_dot = A * PHI;
    PHI_dot = reshape(PHI_dot, [N_phi^2, 1]);


    % diff eq vector
    if(Np == (Nc + Ns))     % no extra parameters
        dx = [zeros(Nc+Ns, 1); ...
            PHI_dot];
    else    % there's extra parameters
        final = Np;
        init = Np - (varObj.Nm*2 - 1);
        % extra parameters time derivative partial
        e_dot = X(init:final).*[0;1;0;1;0;1;0;1;0;1;0;1]; 
        E_dot = [e_dot(2); 0; e_dot(4); 0; e_dot(6); 0; e_dot(8); 0;...
            e_dot(10); 0 ;e_dot(12); 0];
        dx = [zeros(Nc+Ns, 1); ...
              E_dot; PHI_dot];
    end
end

