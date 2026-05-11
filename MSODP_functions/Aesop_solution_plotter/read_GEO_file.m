function [C, S, sigma_C, sigma_S] = read_GEO_file(filename, n_max)
%READ_RECOEF_FILE Read spherical harmonic coefficients from RECOEF file
%
% Assumed ordering in file:
%   (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), ..., (n_max,n_max)
%
% The function DOES NOT use the degree/order written in the file, because
% spacing may be corrupted. Instead, it reads the coefficients in sequence.
%
% From each RECOEF line, it extracts the LAST 5 numeric values:
%   [C, S, sigma_C, sigma_S, extra]
% and ignores the final extra value (typically -1.0).
%
% Outputs are stored as matrices such that:
%   C(n+1,m+1)       = C_nm
%   S(n+1,m+1)       = S_nm
%   sigma_C(n+1,m+1) = sigma_C_nm
%   sigma_S(n+1,m+1) = sigma_S_nm
%
% Entries for undefined coefficients remain NaN.

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    cleaner = onCleanup(@() fclose(fid));

    % Preallocate with NaNs
    C       = NaN(n_max+1, n_max+1);
    S       = NaN(n_max+1, n_max+1);
    sigma_C = NaN(n_max+1, n_max+1);
    sigma_S = NaN(n_max+1, n_max+1);

    % Total number of coefficients from n=2 to n=n_max
    n_coeff = sum((2:n_max) + 1);

    % Regex for integers / decimals / scientific notation
    num_pattern = '[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eEdD][-+]?\d+)?';

    % Current coefficient index in the assumed ordering
    n = 2;
    m = 0;
    count = 0;

    while ~feof(fid) && count < n_coeff
        line = fgetl(fid);

        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end

        % Only process RECOEF lines
        if ~contains(line, 'RECOEF')
            continue;
        end

        % Extract all numeric tokens from the line
        toks = regexp(line, num_pattern, 'match');

        % We expect at least 5 numbers at the end:
        % [C, S, sigma_C, sigma_S, extra]
        if numel(toks) < 5
            warning('Skipping malformed line:\n%s', line);
            continue;
        end

        vals = str2double(toks);

        % Take the last 5 numeric values, ignore the last one
        tail = vals(end-4:end);

        cval   = tail(1);
        sval   = tail(2);
        sigc   = tail(3);
        sigs   = tail(4);

        % Store
        C(n+1, m+1)       = cval;
        S(n+1, m+1)       = sval;
        sigma_C(n+1, m+1) = sigc;
        sigma_S(n+1, m+1) = sigs;

        % Advance (n,m) in the assumed ordering
        count = count + 1;
        m = m + 1;
        if m > n
            n = n + 1;
            m = 0;
        end
    end

    if count < n_coeff
        warning('Only %d of %d coefficients were read up to degree %d.', ...
            count, n_coeff, n_max);
    end
end