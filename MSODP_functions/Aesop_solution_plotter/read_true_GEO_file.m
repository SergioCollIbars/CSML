function [C, S, sigma_C, sigma_S] = read_true_GEO_file(filename, n_max)
%READ_TRUE_RECOEF_FILE Read true SH coefficients from a RECOEF-like file
%
% This reader IGNORES the degree/order written in the file because spacing
% may be corrupted for large values.
%
% It assumes the coefficients appear in this order:
%
%   1) Zonal terms first:
%      (2,0), (3,0), (4,0), ..., (n_max,0)
%
%   2) Then all tesseral/sectorial terms:
%      (2,1), (2,2), (3,1), (3,2), (3,3), (4,1), (4,2), ..., (n_max,n_max)
%
% From each RECOEF line, it extracts the LAST 5 numeric values:
%   [C, S, sigma_C, sigma_S, extra]
% and ignores the final extra value (typically -1.0).
%
% Outputs:
%   C(n+1,m+1)       = C_nm
%   S(n+1,m+1)       = S_nm
%   sigma_C(n+1,m+1) = sigma_C_nm
%   sigma_S(n+1,m+1) = sigma_S_nm
%
% Undefined entries remain NaN.

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    % Preallocate
    C       = NaN(n_max+1, n_max+1);
    S       = NaN(n_max+1, n_max+1);
    sigma_C = NaN(n_max+1, n_max+1);
    sigma_S = NaN(n_max+1, n_max+1);

    % Build the expected (n,m) ordering
    nm_list = zeros(sum((2:n_max)+1), 2);
    k = 1;

    % 1) Zonals first: (2,0), (3,0), ..., (n_max,0)
    for n = 2:n_max
        nm_list(k,:) = [n, 0];
        k = k + 1;
    end

    % 2) Then m>0 terms:
    %    (2,1), (2,2), (3,1), (3,2), (3,3), ...
    for n = 2:n_max
        for m = 1:n
            nm_list(k,:) = [n, m];
            k = k + 1;
        end
    end

    n_coeff = size(nm_list, 1);

    % Regex for numeric values, including scientific notation
    num_pattern = '[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eEdD][-+]?\d+)?';

    count = 0;

    while ~feof(fid) && count < n_coeff
        line = fgetl(fid);

        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end

        if ~contains(line, 'RECOEF')
            continue;
        end

        toks = regexp(line, num_pattern, 'match');

        % Need at least [C, S, sigma_C, sigma_S, extra]
        if numel(toks) < 5
            warning('Skipping malformed line:\n%s', line);
            continue;
        end

        vals = str2double(strrep(toks, 'D', 'E'));
        vals = str2double(strrep(string(toks), 'D', 'E'));
        vals = double(vals);

        tail = vals(end-4:end);

        cval = tail(1);
        sval = tail(2);
        sigc = tail(3);
        sigs = tail(4);

        count = count + 1;
        n = nm_list(count, 1);
        m = nm_list(count, 2);

        C(n+1, m+1)       = cval;
        S(n+1, m+1)       = sval;
        sigma_C(n+1, m+1) = sigc;
        sigma_S(n+1, m+1) = sigs;
    end

    if count < n_coeff
        warning('Only %d of %d coefficients were read up to degree %d.', ...
            count, n_coeff, n_max);
    end
end