function [Y, state] = compute_measurements(instrumentParams, planetParams, ...
    poleParams, time, state, Cnm, Snm, BN_mat)
    
    % extract params
    Nt = length(time);

    GM          = planetParams(1); 
    Re          = planetParams(2);
    n_max       = planetParams(3);
    normalized  = planetParams(4);

    W           = poleParams(1);
    W0          = poleParams(2);
    RA          = poleParams(3);
    DEC         = poleParams(4);

    % generate noise profile
% %     Nm = length(instrumentParams(:, 1));
% %     noise = ones(Nm, Nt) * NaN;
% %     for k = 1:Nm
% %         [n] = generate_gradNoise(instrumentParams(k, :), Nt);
% %         noise(k, :) = n';
% %     end

    % generate white noise
    sigma_m = sqrt(instrumentParams(1, 5)) * instrumentParams(1, 2);
    noise_white = normrnd(0, sigma_m, [6, Nt]);

    % generate flicker (alpha = 2) noise
    S_w        = instrumentParams(1, 2);                           % [mE/sqrt(Hz)]

    fs         =  instrumentParams(1, 5);                          % [Hz]
    f_min      =  instrumentParams(1, 6);                          % [Hz]

    sigma_eta2 = 4 * S_w^2 * sin(pi * f_min / fs)^2;

    eta             = sqrt(sigma_eta2) * randn(6, Nt);    % RW increments
    noise_flicker   = cumsum(eta')';                      % random walk

    %% compute measurements
    Y = ones(6, Nt) * NaN;
    for k = 1:Nt
        Wt = W0 + W * time(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        r = state(k, 1:3)'; v = state(k, 4:6)';
        [Y_ACI, ~] = gradiometer_meas(time(k) ,[GM, Re, n_max, normalized],...
            ACAF_ACI, [r', v'], zeros(9, 1), Cnm, Snm);

        % rotation to Instrument frame
        maxInd = 3 * k; minInd = maxInd - 2;
        BN = BN_mat(minInd:maxInd, :);

        T_ACI = [Y_ACI(1),Y_ACI(2),Y_ACI(3);Y_ACI(4),Y_ACI(5),Y_ACI(6);...
                 Y_ACI(7),Y_ACI(8),Y_ACI(9)];
        T_B   = BN * T_ACI * BN';

        Y(:, k) = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);T_B(2,3);T_B(3,3)]./1E-12;   % [mE]
    end
    
    % add noise
    cnst_bias = instrumentParams(:, 3).*ones(6, Nt);
    Y         = Y + noise_white + cnst_bias + noise_flicker;
    
    % update true bias
    state(:, 7:12) = noise_flicker' + cnst_bias';
end


function [noise] = generate_gradNoise(instrumentParams, Nt)
    % sampling frequency
    fs     = instrumentParams(1, 5);

    % generate noise profile
    k  = (0:Nt-1)';
    f  = k*fs/Nt;

    Swhite = instrumentParams(1, 2);                                        % white noise level [units^2/Hz]
    alpha  = instrumentParams(1, 4);                                        % 1/f^alpha
    fL     = instrumentParams(1, 6);

    % One-sided frequency vector (0..fs/2)
    f1 = f(1:Nt/2+1);          % this is what pwelch uses
    
    S_noise = ones(size(f1)); % base flat level in band and above
    idxLow  = (f1 < fL & f1 > 0);
    
    % 1/f below fL, flat (value = 1) from fL up
    S_noise(idxLow) = (fL ./ (f1(idxLow).^alpha));
    
    % If you want a custom *flat* level, just scale S_noise:
    S_flat = Swhite^2;            % your desired flat PSD value
    S_noise = S_flat * S_noise;
    
    % --------- Build frequency-domain spectrum with correct scaling ---------
    phi = 2*pi*rand(Nt/2+1,1);          % phases for 0..Nyquist
    
    Spec = zeros(Nt,1);
    
    % DC (k=1 in MATLAB)
    Spec(1) = sqrt(S_noise(1) * fs * Nt) .* exp(1j*phi(1));
    
    % Positive freqs (excluding DC and Nyquist)
    kpos = 2:Nt/2;
    Spec(kpos) = sqrt(S_noise(kpos) * fs * Nt / 2) .* exp(1j*phi(kpos));
    
    % Nyquist (k = N/2+1)
    Spec(Nt/2+1) = sqrt(S_noise(end) * fs * Nt) .* exp(1j*phi(end));
    
    % Negative freqs by Hermitian symmetry
    Spec(Nt/2+2:end) = conj(Spec(Nt/2:-1:2));
    
    % Time-domain noise (real)
    noise = ifft(Spec,'symmetric');    % 'symmetric' keeps tiny imag parts away

   
% %     %% PLOT NOISE
% %     nfft_psd = Nt;
% %     [PSD_est, f_est] = pwelch(noise, hamming(nfft_psd/2), nfft_psd/4, ...
% %         nfft_psd, fs, 'onesided');    
% %     
% %     figure;
% %     loglog(f_est, PSD_est, 'LineWidth',1.5); hold on;
% %     loglog(f_est, S_noise, 'LineWidth',1.5); hold on;
% %     xlabel('Frequency [Hz]');
% %     ylabel('PSD [units^2/Hz]');
% %     legend('Estimated PSD (pwelch)','Model PSD');
% %     title('Comparison: model vs estimated PSD');
% %     grid on;
end
