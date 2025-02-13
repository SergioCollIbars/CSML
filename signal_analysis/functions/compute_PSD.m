function [Pxx] = compute_PSD(FFT, Fs, nfft, N)
    % compute_PSD function
    % Description: compute the Power Spectral density (PSD) of the fft
    % frequency series. 
    % Input: FFT: fast fourier transform of the time series x(t)
    %        FS: Frequency sample (Hz)
    %        nfft: number of points of the FFT
    %        N: time series lenght
    % Output: PSD: power spectral density
    Pxx2 = abs(FFT).^2 /Fs /N; 
    Pxx2 = reshape(Pxx2, [nfft, 1]);
    Pxx = [Pxx2(1); 2*Pxx2(2:(nfft/2))];
end