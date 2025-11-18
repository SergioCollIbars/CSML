function [b_GG,b_acc] = load_bias(planetParams)
    % Return the gradiometer bias and accelerometer bias in non-dimensional
    % units.
    % Date: 09/25/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % measurement units
    GG_measDim  = 1E-9;                                  % [s^-2]
    Acc_measDim = 1;                                     % [m/s^2]

    % Initial bias units gradiometer [s^-2]
    Bxx = 0.9E-9;
    Bxy = 4E-9;
    Bxz = 30E-9;
    Byy = 6E-9;
    Byz = 11E-9;
    Bzz = 0.7E-9;
    b_GG = [Bxx;Bxy;Bxz;Byy;Byz;Bzz]./GG_measDim.*1E-3.*0; % [Kilo-Eotvos]

    % Initial accelerometer units gradiometer [m s^-2]
    Bx = 1.5;
    By = .4;
    Bz = 1;
    b_acc = [Bx;By;Bz]./Acc_measDim;                    % [m/s^2]
    
end

