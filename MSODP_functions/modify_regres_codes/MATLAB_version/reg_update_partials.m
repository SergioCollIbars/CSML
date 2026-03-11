function info = reg_update_partials(inReg, outReg, rowIdx, parCols, newVals, opts)
%REG_UPDATE_PARTIALS_ROWS Update selected partial columns for selected observation
% rows in a binary .reg file, copying all other bytes unchanged.
%
% Inputs
%   inReg   : input .reg (binary)
%   outReg  : output .reg (binary)
%   rowIdx  : Mx1 or 1xM observation indices to modify (1-based)
%   parCols : 1xK parameter indices within the partials block (1..npar)
%   newVals : MxK replacement values for those rows/cols
%
% Options (opts)
%   opts.extraDoubles          : number of "extra" doubles before partials (default 5)
%   opts.partialsStartInRecord : 1-based start index of partials inside each record
%                               (default extraDoubles+1)
%   opts.endianness            : 'ieee-le' (default) or 'ieee-be'
%
% Output
%   info struct with npar, rowsize, nObs, etc.

    % ----------------------------
    % Defaults & basic checks
    % ----------------------------
    if nargin < 6 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'extraDoubles') || isempty(opts.extraDoubles), opts.extraDoubles = 5; end
    if ~isfield(opts,'partialsStartInRecord') || isempty(opts.partialsStartInRecord)
        opts.partialsStartInRecord = opts.extraDoubles + 1;
    end
    if ~isfield(opts,'endianness') || isempty(opts.endianness), opts.endianness = 'ieee-le'; end

    assert(isfile(inReg), 'Input file not found: %s', inReg);

    rowIdx = rowIdx(:);
    rowIdx = sort(rowIdx);
    M = numel(rowIdx);
    assert(M >= 1, 'rowIdx is empty.');
    assert(all(diff(rowIdx) > 0), 'rowIdx must be unique (no duplicates).');

    parCols = parCols(:).';
    K = numel(parCols);
    assert(size(newVals,1) == M, 'newVals must have M rows (same as rowIdx).');
    assert(size(newVals,2) == K, 'newVals must have K columns (same as parCols).');

    % ----------------------------
    % Open files
    % ----------------------------
    fin  = fopen(inReg, 'r', opts.endianness);
    assert(fin ~= -1, 'Could not open input file.');

    fout = fopen(outReg, 'w', opts.endianness);
    assert(fout ~= -1, 'Could not open output file.');

    cleanupObj = onCleanup(@() cleanupFiles(fin, fout));

    % ----------------------------
    % Small helpers
    % ----------------------------
    rdChars = @(n) fread(fin, n, 'uint8')';
    wrChars = @(b) fwrite(fout, b, 'uint8');

    rdI64 = @() fread(fin, 1, 'int64');
    wrI64 = @(x) fwrite(fout, x, 'int64');

    rdD = @(n) fread(fin, n, 'double');
    wrD = @(x) fwrite(fout, x, 'double');

    % ----------------------------
    % Copy header (mirrors your C)
    % ----------------------------
    PROGNAM_LNGTH = 10;
    RUNDES_LNGTH  = 280;
    DATIM_LNGTH   = 10;
    STANAM_LNGTH  = 10;
    PARNAM_LNGTH  = 10;

    wrChars(rdChars(PROGNAM_LNGTH));
    wrChars(rdChars(RUNDES_LNGTH));
    wrChars(rdChars(DATIM_LNGTH));
    wrChars(rdChars(DATIM_LNGTH));

    narc = rdI64(); wrI64(narc);
    nsta = rdI64(); wrI64(nsta);

    if nsta ~= 0
        wrChars(rdChars(STANAM_LNGTH * double(nsta)));                % stanam
        idsta = fread(fin, double(nsta), 'int64'); fwrite(fout, idsta, 'int64');
        wrD(rdD(1));                                                  % re
        wrD(rdD(1));                                                  % flat
        wrD(rdD(3 * double(nsta)));                                   % stacor (3 per station)
    end

    ns   = rdI64(); wrI64(ns);
    nf   = rdI64(); wrI64(nf);
    np   = rdI64(); wrI64(np);
    nga  = rdI64(); wrI64(nga);
    iarc = rdI64(); wrI64(iarc);
    wrD(rdD(1));                                                      % utc

    igrps = rdI64(); wrI64(igrps);

    if igrps > 0
        igrpcnt = fread(fin, double(igrps), 'int64'); fwrite(fout, igrpcnt, 'int64');
        igrptot = fread(fin, double(igrps), 'int64'); fwrite(fout, igrptot, 'int64');
    end

    idsat = fread(fin, double(ns), 'int64'); fwrite(fout, idsat, 'int64'); % idsat
    wrD(rdD(6 * double(ns)));                                            % stasat

    npar = 6*double(ns) + double(np);

    wrD(rdD(npar));                                                      % aprval_loc
    wrChars(rdChars(PARNAM_LNGTH * npar));                               % parnam_loc (10 chars each)
    wrD(rdD(npar));                                                      % parsig_loc
    wrD(rdD(npar));                                                      % parsca_loc

    % ----------------------------
    % Observation records
    % ----------------------------
    extra   = double(opts.extraDoubles);
    rowsize = npar + extra + double(igrps);

    assert(all(parCols >= 1 & parCols <= npar), 'parCols must be within 1..npar (npar=%d).', npar);

    p0 = double(opts.partialsStartInRecord); % 1-based index within record
    assert(p0 >= 1 && (p0 + npar - 1) <= rowsize, ...
        'partialsStartInRecord=%d invalid for rowsize=%d and npar=%d.', p0, rowsize, npar);

    partialBlock = p0 : (p0 + npar - 1);
    overwriteIdx = partialBlock(parCols);

    obs = 0;
    targetPtr = 1;

    while true
        rec = fread(fin, rowsize, 'double');
        if numel(rec) < rowsize
            break; % EOF
        end
        obs = obs + 1;

        if targetPtr <= M && obs == rowIdx(targetPtr)
            rec(overwriteIdx) = newVals(targetPtr,:).';
            targetPtr = targetPtr + 1;
        end

        fwrite(fout, rec, 'double');
    end

    assert(targetPtr == M+1, ...
        'File ended before reaching requested rowIdx=%d (processed obs=%d).', rowIdx(targetPtr), obs);

    info = struct('npar',npar, 'rowsize',rowsize, 'igrps',double(igrps), ...
                  'extra',extra, 'partialsStartInRecord',p0, ...
                  'nObs',obs, 'modifiedRows',rowIdx, 'parCols',parCols);
end

function cleanupFiles(fin, fout)
    if ~isempty(fin)  && fin  ~= -1, fclose(fin);  end
    if ~isempty(fout) && fout ~= -1, fclose(fout); end
end