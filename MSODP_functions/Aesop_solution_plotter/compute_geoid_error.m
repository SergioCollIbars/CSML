function [degree_error_mm,degree_sigma_mm] = compute_geoid_error(Nmax,...
            C_est, S_est, C_ref, S_ref, sigC, sigS, R)
    % Inputs:
    % C_est, S_est : estimated normalized coefficients
    % C_ref, S_ref : reference normalized coefficients
    % sigC, sigS   : coefficient formal uncertainties
    % R            : reference radius [m]
    % Nmax         : maximum degree
    
    degree_error_mm = zeros(Nmax+1,1);
    degree_sigma_mm = zeros(Nmax+1,1);
    
    for n = 0:Nmax
        err_sum = 0;
        sig_sum = 0;
    
        for m = 0:n
            dC = C_est(n+1,m+1) - C_ref(n+1,m+1);
            dS = S_est(n+1,m+1) - S_ref(n+1,m+1);
    
            err_sum = err_sum + dC^2 + dS^2;
            sig_sum = sig_sum + sigC(n+1,m+1)^2 + sigS(n+1,m+1)^2;
        end
    
        degree_error_mm(n+1) = 1000 * R * sqrt(err_sum);
        degree_sigma_mm(n+1) = 1000 * R * sqrt(sig_sum);
    end
end