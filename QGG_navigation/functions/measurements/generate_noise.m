function [noise] = generate_noise(sigma, Nt, Nn, flicker, Fs)
    % output value    
    noise = ones(length(sigma), Nt);

    % generate flicker noise or gaussian
    if(flicker)
        for j = 1:length(sigma)
            measDim = 7.1040e-12;       % [1/s^2]

            sigma_wn = sigma(j) * measDim;
            s =  6.32E-12 / sigma_wn;   % scaling value
            sigma_pn = sigma_wn / (50 * s);
        
            % generate noise [1/s^2]
            [X, ~] = noise_profile(sigma_wn, sigma_pn, Nt, Fs);
        
            % save noise and PSD
            X = X./measDim; % [-]
            noise(j, :)  = X * Nn;
        end
    else
        for j=1:length(sigma)
            noise(j, :) = normrnd(0, sigma(j), [1, Nt]) * Nn;
        end
    end
end

