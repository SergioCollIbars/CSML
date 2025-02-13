function [state_true, state_recons, state_error, TIME] =...
    orbit_corrector...
    (X0, planetParams, poleParams ,Cmat, Smat, tt, enviroment, order)
    

    % define integration options
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    STM_0 = reshape(eye(6,6), [36, 1]);

    % ODE 113
    [~, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat, Smat, enviroment), tt, [X0;STM_0], options);
    t = tt;
    TIME = tt;
    
    % dimensionalize equations and rotate to inertial frame
    if(enviroment == "CR3BP")
        ddU = ones(9, length(t)) * NaN;
        for j = 1:length(t) % generating measurements
            m_1 = 5.974E24;  % [Kg]
            m_2 = 7.348E22;  % [Kg]
            pi1 = m_1 / (m_1 + m_2);
            pi2 = m_2 / (m_1 + m_2);

            r_B = state_true(j, 1:3)';

            r1_B = r_B + [+pi2; 0; 0];
            r2_B = r_B + [-pi1; 0; 0];

            [T] = gradmeas_rotFrame(pi2, r_B(1), r_B(2), r_B(3), ...
                vecnorm(r1_B), vecnorm(r2_B));

            ddU(:, j) = reshape(T, [9, 1]); % rotating frame. Non-dim units
        end
    elseif(enviroment == "2BP")
        % generate gradiometer measurements
        [~, ddU] = gradiometer_meas(t, state_true(:, 1:3), Cmat, Smat, ...
            poleParams, planetParams);

    elseif(enviroment == "3BP")
        ddU = ones(9, length(t)) * NaN;
        J = zeros(1, length(t));
        acc = zeros(3, length(t));
        
        % planets params
        mu = planetParams(1);       % [-]

        % rotate to inertial frame
        [r_I, v_I] = rotate2inertial(state_true(:, 1:3)', state_true(:, 4:6)',...
            t, 1);

        % Jacobi integral
        r_B = state_true(:, 1:3)';
        v_B = state_true(:, 4:6)';
        r1 = r_B - [-mu;0;0];
        r2 = r_B - [1-mu;0;0];

        Jb = -0.5.*(vecnorm(r_B).^2) - ...
            ((1-mu)./vecnorm(r1) + mu./vecnorm(r2)) + 0.5.*(vecnorm(v_B).^2);
        figure()
        plot(t, Jb)

        % compute gradiometer measurements
        for j = 1:length(t)
            % Earth and Moon motion (position and velocity)
            re_B = [-mu;0;0];
            ve_B = [0;-mu;0];
            rm_B = [(1-mu);0;0];
            vm_B = [0; (1-mu);0];

            [re_I, ~] = rotate2inertial(re_B, ve_B, t(j), 1);
            [rm_I, ~] = rotate2inertial(rm_B, vm_B, t(j), 1);

            % relative position to the spacecraft
            res = r_I(:, j) - re_I;  % [m]
            rms = r_I(:, j) - rm_I;  % [m]

            % compute SOGT
            T = gradmeas_CR3BP_inertial(mu,r_I(1, j), r_I(2, j),...
                r_I(3, j), t(j));
            
            % compute gravity acc
            acc(:, j) = - mu/(vecnorm(rms)^3) * rms - (1-mu) /(vecnorm(res)^3) * res;
            
            % Energy in the inertial frame
            U = 1*(mu/vecnorm(rms) + (1-mu)/vecnorm(res)); 

            H = cross(r_I(:, j), v_I(:, j));
            W = [0;0;1];

            J(j) = 0.5 * vecnorm(v_I(:, j))^2 - dot(W, H) - U;

            ddU(:, j) = reshape(T, [9, 1]); % no rotating frame. Non-dim units
        end
        figure()
        plot(t, J)

        state_true(:, 1:3) = r_I';
        state_true(:, 4:6) = v_I';
    end
    
    % reconstruct trajectory
    delta = (state_true(2, 1:6) - state_true(1, 1:6))';
    [delta_new] = reconstruct_traj(t, ddU, delta, enviroment, order);
    
    % state reconstruction from nominal
    state_recons = zeros(6, length(t));
    state_recons(:, 1) = state_true(1,1:6)';
    for j = 2:length(t)
% %         PHI_new = reshape(state_true(j, 7:end), [6,6]);
% %         delta_new2  = PHI_new * delta;
% %         state_recons(:, j) = state_recons(:, j-1) + delta_new2;
        state_recons(:, j) = state_recons(:, j-1) + delta_new(:, j);
    end
    
    % error states
    state_error = state_true(:, 1:6) - state_recons';
end