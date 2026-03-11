function [x] = create_noise_from_PSD(fs, T, plotCommand)
    
    % degree = 0; % integration PSD degree (n) (n>0 : derivative) (n<0: integral)
    %% Sampling parameters
    N  = fs*T;         % samples (make sure N is even)
    assert(mod(N,2)==0, 'N must be even.');
    df = fs/N;
    fpos = (0:N/2).' * df;     % one-sided freq grid: 0 ... fs/2

    %% Target one-sided PSD S1(f)  [units^2/Hz]
    % Example: S1(f) = sqrt(f)
    % S0 = 8.5E-6 * sqrt(fpos.^(-1));
    
    S1        = 3E-8 * sqrt(1 + 4.6E-8 * fpos.^(-2));
    S2_ASTRIX = (((2*pi*fpos).^(-1)).*S1).^2;

    S2_ST     = (8.5E-6 * sqrt(fpos.^(-1))).^2;
    
    S0        = 1E-10 * sqrt(0.4 + 0.001 * fpos.^(-1) + 2500*fpos.^(4));
    S2_ACC    = (((2*pi*fpos).^(-2)).*S0).^2;

    S1 = (S2_ACC.^(-1) + S2_ST.^(-1) + S2_ASTRIX.^(-1)).^(-1);
    S1(1) = 0;                % set DC PSD to 0 for this example (avoid DC blow-up)
    
    % (Optional) set Nyquist PSD too:
    S1(end) = S1(end);
    
    %% Build FFT bins X such that E{|X(k)|^2} = fs*N*S2(fk)
    % where S2 = S1/2 for positive freqs (excluding DC and Nyquist)
    %rng(1);
    
    X = zeros(N,1);
    
    % DC bin (k=1) must be real
    X(1) = sqrt(fs*N*S1(1)) * randn;
    
    % Positive freqs k = 2 ... N/2 (exclude Nyquist)
    kpos = (2:N/2).';
    fk   = fpos(kpos);
    
    % Complex Gaussian with unit variance: (a+jb)/sqrt(2)
    G = (randn(length(kpos),1) + 1i*randn(length(kpos),1)) / sqrt(2);
    
    % Scale to match PSD:
    % For real signal, two-sided S2 = S1/2 at these bins
    X(kpos) = G .* sqrt(fs*N*(S1(kpos)/2));
    
    % Nyquist bin (k=N/2+1) must be real
    knyq = N/2 + 1;
    X(knyq) = sqrt(fs*N*S1(end)) * randn;
    
    % Negative freqs: Hermitian symmetry
    X(N/2+2:end) = conj(flipud(X(2:N/2)));
    
    %% Time-domain realization
    x = real(ifft(X));         % MATLAB ifft includes 1/N internally
    t = (0:N-1)'/fs;
    
    if(plotCommand)
        figure;
        plot(t, x);
        xlabel('Time [s]'); ylabel('x(t)');
        title('Time-domain noise realization');
        grid on;
    end
    
    %% Verify PSD using Welch (one-sided for real x)
    nfft = 4096;
    [pxx, fwelch] = pwelch(x, hanning(2048), 1024, nfft, fs, 'onesided');
    
    % Interpolate target PSD onto Welch frequencies
    S1_w = interp1(fpos, S1, fwelch, 'linear', 'extrap');
    
    if(plotCommand)
        figure;
        loglog(fwelch, pxx, 'LineWidth', 1.5); hold on;
        loglog(fwelch, S1_w, '--', 'LineWidth', 1.5);
        xlabel('Frequency [Hz]'); ylabel('PSD [units^2/Hz]');
        legend('Estimated PSD (Welch)','Target PSD');
        title('PSD Comparison');
        grid on;
    end
    
    %% Optional sanity check: variance consistency
    % For a one-sided PSD, variance ~= integral_0^{fs/2} S1(f) df
    var_from_psd = trapz(fpos, S1);
    var_time     = var(x);
    
    fprintf('Var from PSD integral: %.6e\n', var_from_psd);
    fprintf('Var from time series : %.6e\n', var_time);
end

