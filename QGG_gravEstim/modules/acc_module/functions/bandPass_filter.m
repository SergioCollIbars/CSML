function [BP_FFT] = bandPass_filter(FFT, freq, MB)
    % Decription: Bandpass the fourier signal for the values specified in 
    %  the MB vector.
    % Input: FFT: values of th Fourier Transform
    %        freq: frequency vector
    %        MB: 2x1 vector with the band passs values
    % Ouput: BP_FFT: values of the band pass Fourier Transform

    ind1 = find((MB(1) < freq) & (freq < MB(2)));
    ind2 = length(FFT) + 1 - ind1;
    
    % Keep values inside the MB
    BP_FFT = zeros(length(FFT), 1);
    BP_FFT(ind1) = FFT(ind1);
    BP_FFT(ind2) = FFT(ind2);
end
