function [Cnm, Snm, sigmaC, sigmaS, n_max] = read_GRACE_GSM(filename, n_max)

    % READ_GRACE_GSM reads GRACE/GRACE-FO GSM spherical harmonic files.
    %
    % Outputs:
    %   Cnm    : (n_max+1) x (n_max+1) cosine coefficients
    %   Snm    : (n_max+1) x (n_max+1) sine coefficients
    %   sigmaC : (n_max+1) x (n_max+1) formal uncertainty of Cnm
    %   sigmaS : (n_max+1) x (n_max+1) formal uncertainty of Snm
    %
    % Indexing:
    %   Cnm(n+1,m+1) corresponds to degree n, order m.

    fid = fopen(filename, 'r');

    if fid < 0
        error('Could not open file: %s', filename);
    end

    % If n_max is not provided, read it from the header
    if nargin < 2 || isempty(n_max)
        n_max = [];

        while ~feof(fid)
            line = fgetl(fid);

            if contains(line, 'degree')
                tokens = regexp(line, 'degree\s*:\s*(\d+)', 'tokens');

                if ~isempty(tokens)
                    n_max = str2double(tokens{1}{1});
                    break
                end
            end
        end

        if isempty(n_max)
            fclose(fid);
            error('Could not determine n_max from file header.');
        end

        frewind(fid);
    end

    Cnm    = zeros(n_max+1, n_max+1);
    Snm    = zeros(n_max+1, n_max+1);
    sigmaC = zeros(n_max+1, n_max+1);
    sigmaS = zeros(n_max+1, n_max+1);

    while ~feof(fid)
        line = fgetl(fid);

        if startsWith(strtrim(line), 'GRCOF2')
            data = textscan(line, '%s %d %d %f %f %f %f %s %s %s');

            n = data{2};
            m = data{3};

            if n <= n_max && m <= n_max
                Cnm(n+1,m+1)    = data{4};
                Snm(n+1,m+1)    = data{5};
                sigmaC(n+1,m+1) = data{6};
                sigmaS(n+1,m+1) = data{7};
            end
        end
    end

    fclose(fid);
end