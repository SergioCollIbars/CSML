function [X ,PSD, frec] = noise_profile(sigma_wn, sigma_pn, npts, Fs)
    % Description: compute the noise profile using white and pink noise
    % inputs: sigma_wn: STD white noise
    %         sigma_pn: STD pink noise
    %         npts: number of points
    %         Fs: sampling frequency
    % outputs: X: noise time series
    %          PSD: power spectral density
    %          frec: frequency series
    
    % number of points FFT
    nfft = 2^nextpow2(npts);
    frec = 0 : Fs/(nfft-1) : Fs/2;

    % Gaussian noise time series
    X_wn = sigma_wn.*randn(npts,1);

    % Pink noise time series. FIR method
    [X_pn] = noiseGenerator_FIR(npts, sigma_pn, 1);
    
    % Band pass frequencies
    MB = [5E-3, 0.1];
    %MB = [5E-5, 0.1];
    [BP_FFT_wn] = bandPass_filter(fft(X_wn, npts), frec, MB);
    MB = [0, 5E-3];
    %MB = [0 , 5E-5];
    [BP_FFT_pn] = bandPass_filter(fft(X_pn, npts), frec, MB);
    BP_FFT = BP_FFT_pn + BP_FFT_wn;
    
    % compute noise from FFT to time series
    X = ifft(BP_FFT, npts);
    X = real(X);

    % total noise PSD
    PSD = compute_PSD(fft(X, nfft), Fs, nfft, npts);   
end

