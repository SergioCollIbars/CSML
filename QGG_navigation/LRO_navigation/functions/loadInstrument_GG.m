function [instrumentParams] = loadInstrument_GG(folder_Name)
    % Description: give the Gravity Gradiometer (GG) instrument parameters
    p = readParams("data/"+folder_Name+"/Instrument_GG.txt");
    
    % Measurement mask
    Mask = [p.mask_xx;p.mask_xy;p.mask_xz;...
            p.mask_yy;p.mask_yz;p.mask_zz];      

    % measurement noise: [mE * sqrt(Hz)]
    sigma = [p.sigma_xx;p.sigma_xy;p.sigma_xz;...
             p.sigma_yy;p.sigma_yz;p.sigma_zz];         

    % Measurement band [Hz]
    fmin = ones(6, 1).*p.fmin;
    fmax = ones(6, 1).*p.fmax;

    % noise slope
    alpha = ones(6, 1).*p.alpha;

    % measurement bias [mE]
    bias = [p.bias_xx;p.bias_xy;p.bias_xz;...
            p.bias_yy;p.bias_yz;p.bias_zz];

    % Scale Factor (SF) [-]
    SF = [p.SF_xx;p.SF_xy;p.SF_xz;...
            p.SF_yy;p.SF_yz;p.SF_zz];

    % sampling frequency [Hz]
    fs = ones(6,1).*p.fs;

    instrumentParams = [Mask, sigma, bias, alpha, fs, fmin, fmax, SF];
end

