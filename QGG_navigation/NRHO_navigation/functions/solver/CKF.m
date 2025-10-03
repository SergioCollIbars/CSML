function [X_hat, P, P_bar, Qn] = CKF(dY, Hi, R0, P, X_hat, PHI, Q)
    %%                     CLASSIC KALMAN FILTER FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 02/06/2023
    %
    %   Description: Compute normal equation at each time state.
    %   using Kalman Filter
    %
    %   Input: 
    %       rho: data range measurements
    %       rhoDot: data range rate measurement
    %       rho_c: computed range measurement
    %       rhoDot_c: computed range rate measurement
    %       H: Visibility matrix related to t: i
    %       R0: measurements covariance
    %       P: state covariance at t: i-1
    %       X_hat: state deviation at t: i-1
    %       PHI: STM at t: i-1 / t: i
    %
    %   Output:
    %       X_hat: state deviation at t: i
    %       P: state covariance at t: i
    %       P_bar: state covariance time update at t: i
    % --------------------------------------------------------------------%
  
    
    % values time update
    X_bar = PHI * X_hat;
    P_bar = PHI * P * PHI' + Q;

    % stations number.
    Nx = length(X_hat);

    if(isnan(dY))
        % update meas
        X_hat = X_bar;
        P = P_bar;
        Qn = Q;
    else
        % compute gain
        K = P_bar * Hi'/(Hi * P_bar * Hi' + R0);
        
        % update meas
        X_hat = X_bar + K * (dY - Hi * X_bar);
        P = (eye(Nx, Nx) - K * Hi) * P_bar * (eye(Nx, Nx) - K * Hi)' + ...
            K * R0 * K';
        
        alpha  = 0;
        Qn = alpha * Q + (1-alpha)*K*dY*(dY'*K');
    end
end