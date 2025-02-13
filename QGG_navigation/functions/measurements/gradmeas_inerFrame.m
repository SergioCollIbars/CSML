function [dU_ACI, ddU_ACI] = gradmeas_inerFrame(TIME, R_ACI, C_mat, S_mat, ...
    poleParams, planetParams)
%%                         GRADIOEMTER MEASUREMENTS                      %%
%                                                                         %   
%   Author: Sergio Coll Ibars                                             %
%   Date: 01/31/2024                                                      %
%                                                                         %
%   Description: function to generate gradiometer measurements at each    %
%  time instant.                                                          %
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
    ddU_ACI = ones(9, Nt) * NaN;
    dU_ACI = ones(3, Nt) * NaN;
    
    % Get orbit data
    GM = planetParams(1);
    Re = planetParams(2);
    n_max = planetParams(3);
    Normalized = planetParams(4);

    W = poleParams(1);
    W0 = poleParams(2);
    RA = poleParams(3);
    DEC = poleParams(4);
    
    % make sure it is stack vertically
    R_ACI = reshape(R_ACI, [Nt, 3]);
    
    % sigma error
    sigma = 0; % [1/s^2]

    for j = 1:Nt
        % ACI position vector
        r_ACI = [R_ACI(j, 1); R_ACI(j, 2); R_ACI(j, 3)];

        % Inertial to ACAF rotation 
        Wt = W * TIME(j) + W0;
        ACAF_N = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % ACAF position vector
        r_ACAF = ACAF_N * r_ACI;
        
        % compute S/C acc
        [~, a, T] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                    r_ACAF, Re, GM, ...
                                                    Normalized);
        % rotate from ECEF 2 ECI
        T = ACAF_N' * T * ACAF_N;
        a = ACAF_N' * a;

        % noise generator. Gaussian
        r = normrnd(0, sigma, [9, 1]);
        
        % stack measurements
        ddU_ACI(:, j) = reshape(T, [9, 1]) + r; 
        dU_ACI(:, j) = a;
    end
end

