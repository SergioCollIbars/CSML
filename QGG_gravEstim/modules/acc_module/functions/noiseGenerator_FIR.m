function [X] = noiseGenerator_FIR(npts, sigma2, alpha)
    % noise generator function
    % Description: creates noise with PSD as 1/f^alpha. Algorithm taken
    % from Kasdin 1995
    % Inputs: npts: number of points of the FFT
    %         sigma2: variance of the white noise
    %         alpha: exponential value of 1/f^alpha
    nn = 2*npts;
    ha = alpha/2;
    Qd = sqrt(sigma2);

    hfa = zeros(1, nn);
    wfa = zeros(1, nn);
    hfa(1) = 1;
    wfa(1) = Qd * randn();

    for i = 2:npts
        hfa(i) = hfa(i-1) * (ha+(i-2))/(i-1);
        wfa(i) = Qd * randn();
    end
    % the rest of the sequence is 0 by definition
    
    % perform FFT
    hfa = fft(hfa, npts);
    wfa = fft(wfa, npts);

    % multiply both
    wfa = wfa.*hfa;

    % do inverse FTT
    wfa = ifft(wfa, npts);
    X = wfa./npts;
end
