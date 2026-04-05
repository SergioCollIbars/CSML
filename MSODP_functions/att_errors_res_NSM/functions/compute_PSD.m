function [f, PSD] = compute_PSD(x, dt)

    % x  : 1xN or Nx1 vector (signal in radians)
    % dt : sampling time (seconds)

    fs = 1/dt;                 % sampling frequency (Hz)

    x = x(:);                  % ensure column vector
    x = x - mean(x);           % remove DC component

    % Welch PSD estimate
    nfft = 2^nextpow2(length(x));
    window = hanning(floor(length(x)/8));
    noverlap = floor(length(window)/2);

    [PSD,f] = pwelch(x, window, noverlap, nfft, fs);
end