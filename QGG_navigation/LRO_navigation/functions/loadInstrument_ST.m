function [instrumentParams] = loadInstrument_ST(folder_Name)
    % Description: give the Star Tracket (ST) instrument parameters
    p = readParams("data/"+folder_Name+"/Instrument_ST.txt");   

    % sampling frequency [Hz]
    fs = p.fs;

    % measurement noise: [rad]
    sigma = p.sigma / sqrt(p.fs) * pi / (180 * 3600);         

    % measurement bias [rad]
    bias = p.bias * pi / (180 * 3600);

    instrumentParams = [sigma, bias, fs];
end

