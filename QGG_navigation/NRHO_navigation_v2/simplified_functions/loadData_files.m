function [tmin, X0_true, X0_filter, P0, gamma] = loadData_files()
    % Load the initital conditions from the saved file and generate the
    % initial filter conditions with the given covariance.
    % Date: 10/01/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load files with initial conditions
    path = "GG_navigation_final_params.mat";
    load(path);
    
    path = "gamma.mat";
    load(path);
    Ns = length(X_estim_final);

    % assign parameters
    tmin    = t_final;
    P0      = reshape(P_final, [Ns, Ns]);
    X0_true = X_true_final(1:Ns)';
    
    % generate error given the new covariance
    P0 = round(P0 * 1e24) / 1e24;
    C = (P0 + P0.')/2;
    error = mvnrnd(zeros(Ns, 1), C, 1)';

     X0_filter = X0_true + error;
% %     X0_filter = X_estim_final;
end

