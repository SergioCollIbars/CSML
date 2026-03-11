function H = read_reg_files_partials(regFile, endianness)
%READ_REG_ALL
% Returns a matrix:
%   H = [partials  residual  sigma]
%
% Size:
%   H = Nobs x (npar + 2)

    if nargin < 2
        endianness = 'ieee-le';  % try 'ieee-be' if numbers look wrong
    end

    fin = fopen(regFile, 'r', endianness);
    assert(fin ~= -1, 'Could not open file.');

    % ----------------------------
    % Constants (from your C code)
    % ----------------------------
    PROGNAM_LNGTH = 10;
    RUNDES_LNGTH  = 280;
    DATIM_LNGTH   = 10;
    STANAM_LNGTH  = 10;
    PARNAM_LNGTH  = 10;

    extra = 5;

    % ----------------------------
    % Skip file header
    % ----------------------------
    fread(fin, PROGNAM_LNGTH, 'uint8');
    fread(fin, RUNDES_LNGTH,  'uint8');
    fread(fin, DATIM_LNGTH,   'uint8');
    fread(fin, DATIM_LNGTH,   'uint8');

    fread(fin, 1, 'int64');  % narc
    nsta = fread(fin, 1, 'int64');

    if nsta ~= 0
        fread(fin, STANAM_LNGTH * double(nsta), 'uint8');
        fread(fin, double(nsta), 'int64');
        fread(fin, 1, 'double');
        fread(fin, 1, 'double');
        fread(fin, 3 * double(nsta), 'double');
    end

    % ----------------------------
    % Arc header
    % ----------------------------
    ns  = fread(fin, 1, 'int64');
    fread(fin, 1, 'int64');  % nf
    np  = fread(fin, 1, 'int64');
    fread(fin, 1, 'int64');  % nga
    fread(fin, 1, 'int64');  % iarc
    fread(fin, 1, 'double'); % utc

    igrps = fread(fin, 1, 'int64');

    if igrps > 0
        fread(fin, double(igrps), 'int64');
        fread(fin, double(igrps), 'int64');
    end

    fread(fin, double(ns), 'int64');
    fread(fin, 6 * double(ns), 'double');

    % ----------------------------
    % Parameter definitions
    % ----------------------------
    npar = 6*double(ns) + double(np);

    fread(fin, npar, 'double');
    fread(fin, PARNAM_LNGTH*npar, 'uint8');
    fread(fin, npar, 'double');
    fread(fin, npar, 'double');

    % ----------------------------
    % Observation records
    % ----------------------------
    rowsize = npar + extra + double(igrps);

    raw = fread(fin, Inf, 'double');
    fclose(fin);

    nObs = floor(length(raw) / rowsize);
    raw  = raw(1:nObs*rowsize);

    rec = reshape(raw, rowsize, nObs).';

    % ----------------------------
    % Extract data
    % ----------------------------

    partials  = rec(:, extra+1 : extra+npar);
    residuals = rec(:, 4);
    sigma     = rec(:, 5);

    % Combine into one matrix
    H = [partials residuals sigma];
end