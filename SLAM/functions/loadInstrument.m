function [instrumentParams] = loadInstrument()
    % Description: give the instrument parameters
    
    % Measurement mask
    Mask = [1;1;1;1;1;1];           % xx, xy, xz, yy, yz, zz

    % measurement noise
    sigma = ones(6, 1).*1;          % [mE * sqrt(Hz)]

    % Measurement band [Hz]
    fmin = ones(6, 1).*1E-4;
    fmax = ones(6, 1).*1;

    % noise slope
    alpha = ones(6, 1);

    % measurement bias [mE]
    bias = ones(6, 1).*(0.1);

    % sampling frequency [Hz]
    fs = ones(6,1).*(1/10);

    instrumentParams = [Mask, sigma, bias, alpha, fs, fmin, fmax];
end

