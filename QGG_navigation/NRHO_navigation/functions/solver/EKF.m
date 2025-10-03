function [X_hat, P, Qn] = EKF(dY, Hi, R0, P, PHI, Q)
    %%               EXTENDED KALMAN FILTER FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 02/08/2023
    %
    %   Description: Compute normal equation at each time state.
    %   using Extended Kalman Filter
    %
    %   Input: 
    %       rho: data range measurements
    %       rhoDot: data range rate measurement
    %       rho_c: computed range measurement
    %       rhoDot_c: computed range rate measurement
    %       H: Visibility matrix related to t: i
    %       R0: measurements covariance
    %       P: state covariance at t: i-1
    %       PHI: STM at t: i-1 / t: i
    %
    %   Output:
    %       X_hat: state deviation at t: i
    %       P: state covariance at t: i
    %       dY: prefit vector
    %       Hi: visibility matrix related to t: i
    % --------------------------------------------------------------------%
    
    % values time update
    Nx = length(diag(P));

    X_bar = zeros(Nx, 1);
    P_bar = PHI * P * PHI' + Q;

    if(isnan(dY))
        % update meas
        X_hat = X_bar;
        P = P_bar;
        Qn = Q.*0;
    else
        % compute gain
        K = P_bar * Hi'/(Hi * P_bar * Hi' + R0);

        % update meas
        X_hat = K * dY;
        P = (eye(Nx, Nx) - K * Hi) * P_bar * (eye(Nx, Nx) - K * Hi)' + ...
            K * R0 * K';

        alpha  = 0.5;
        Qn = alpha * Q + (1-alpha)*K*dY*(dY'*K');
    end

end

