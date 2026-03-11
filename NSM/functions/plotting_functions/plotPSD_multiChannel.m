function plotPSD_multiChannel(X, Fs, pltt, lgtt, window, noverlap, nfft)
% plotPSD_multiChannel  Plot PSD for each channel (row) of X in one figure.
%
% INPUTS
%   X        : N x Nt data matrix (N channels, Nt time samples)
%   Fs       : sampling frequency [Hz]
%   window   : (optional) pwelch window (scalar length or vector). Default: Hann(1024)
%   noverlap : (optional) number of overlapping samples. Default: 50% of window length
%   nfft     : (optional) FFT length. Default: max(256, 2^nextpow2(window length))
%
% NOTES
%   - Uses Welch PSD estimate (pwelch)
%   - Plots all channels on a single axes (log-log)
%   - Assumes each row is a time series

    narginchk(2,5);

    if isempty(X) || ~ismatrix(X)
        error('X must be a non-empty 2D matrix of size N x Nt.');
    end
    if ~isscalar(Fs) || Fs <= 0
        error('Fs must be a positive scalar (Hz).');
    end

    [N, Nt] = size(X);

    % Defaults
    if nargin < 5 || isempty(window)
        winLen = min(1024, Nt);
        window = hann(winLen, 'periodic');
    else
        if isscalar(window)
            winLen = min(window, Nt);
            window = hann(winLen, 'periodic');
        else
            winLen = min(length(window), Nt);
            window = window(:);
            if length(window) ~= winLen
                window = window(1:winLen);
            end
        end
    end

    if nargin < 6 || isempty(noverlap)
        noverlap = floor(0.5 * winLen);
    else
        noverlap = min(noverlap, winLen-1);
    end

    if nargin < 7 || isempty(nfft)
        nfft = max(256, 2^nextpow2(winLen));
    end

    % Compute PSDs
    % Store first PSD to initialize plot with correct f vector
    figure('Color','w'); clf;
    ax = axes; hold(ax, 'on'); grid(ax, 'on');

    for i = 1:N
        xi = X(i,:);
        if all(~isfinite(xi))
            warning('Channel %d contains no finite values; skipping.', i);
            continue;
        end

        % Remove mean (helps with DC dominating)
        xi = xi - mean(xi, 'omitnan');

        % Replace NaNs/Infs (pwelch can't handle NaNs)
        bad = ~isfinite(xi);
        if any(bad)
            % simple linear interpolation over bad samples
            t = 1:Nt;
            good = ~bad;
            if nnz(good) < 2
                warning('Channel %d has too few valid samples; skipping.', i);
                continue;
            end
            xi(bad) = interp1(t(good), xi(good), t(bad), 'linear', 'extrap');
        end

        [Pxx, f] = pwelch(xi, window, noverlap, nfft, Fs, 'onesided');

        % Plot (use log-log when possible)
        lw = 1.5;
        if all(f > 0) && all(Pxx > 0)
            loglog(ax, f, sqrt(Pxx), 'DisplayName', lgtt(i),...
                'LineWidth', lw);
        else
            % fallback
            plot(ax, f, sqrt(Pxx), 'DisplayName', lgtt(i), ...
                'LineWidth', lw);
            set(ax, 'YScale', 'log');
            set(ax, 'XScale', 'log');
        end
    end

    xlabel(ax, 'Frequency [Hz]');
    ylabel(ax, '$PSD^{1/2}$ $[E/\sqrt{Hz}]$', 'Interpreter','latex');
    title(ax, pltt);
    legend(ax, 'show', 'Location', 'bestoutside');
end