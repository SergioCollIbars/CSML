function info = reg_update_partials_allRows_v2(inReg, outReg, parCols, newVals, opts)
%REG_UPDATE_PARTIALS_ALLROWS Overwrite selected partial columns AND (residual,sigma)
% in a binary .reg file, preserving everything else.
%
% info = reg_update_partials_allRows(inReg, outReg, parCols, newVals)
%
% Inputs
%   inReg   : input .reg (binary)
%   outReg  : output .reg (binary)
%   parCols : 1xK parameter indices (1..npar) to overwrite in each obs record
%   newVals : Nobs x (K+2) matrix:
%             columns 1:K   -> values to write into partials(:, parCols)
%             column  K+1   -> residual  (writes to rec(4))
%             column  K+2   -> sigma     (writes to rec(5))
%
% opts (optional)
%   opts.extraDoubles          : default 5 (matches your C code)
%   opts.partialsStartInRecord : default extraDoubles+1 (record=[extra][partials][igrps])
%   opts.endianness            : default 'ieee-le'
%   opts.formatCheckOnly       : if true, reads header and stops (no write)

    arguments
        inReg (1,:) char
        outReg (1,:) char
        parCols (1,:) double {mustBeInteger, mustBePositive}
        newVals double
        opts.extraDoubles (1,1) double {mustBeInteger, mustBeNonnegative} = 5
        opts.partialsStartInRecord (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        opts.igrpsAtEnd (1,1) logical = true 
        opts.endianness (1,:) char = 'ieee-le'
        opts.formatCheckOnly (1,1) logical = false
    end

    % Dependent default
    if opts.partialsStartInRecord == 0
        opts.partialsStartInRecord = opts.extraDoubles + 1;
    end

    assert(isfile(inReg), "Input file not found: %s", inReg);

    fin  = fopen(inReg, 'r', opts.endianness);
    assert(fin~=-1, "Could not open input file.");

    if ~opts.formatCheckOnly
        fout = fopen(outReg, 'w', opts.endianness);
        assert(fout~=-1, "Could not open output file.");
    else
        fout = -1;
    end

    cleanupObj = onCleanup(@() cleanupFiles(fin, fout));

    % Helpers
    rdChars = @(n) fread(fin, n, '*uint8')';
    wrChars = @(bytes) fwriteIf(fout, bytes, 'uint8');

    rdI64 = @() fread(fin, 1, '*int64');
    wrI64 = @(x) fwriteIf(fout, x, 'int64');

    rdD  = @(n) fread(fin, n, '*double');
    wrD  = @(x) fwriteIf(fout, x, 'double');

    % -------------------------
    % 1) Copy File Header
    % -------------------------
    PROGNAM_LNGTH = 10;
    RUNDES_LNGTH  = 280;
    DATIM_LNGTH   = 10;
    STANAM_LNGTH  = 10;
    PARNAM_LNGTH  = 10;

    prognam = rdChars(PROGNAM_LNGTH); wrChars(prognam);
    rundes  = rdChars(RUNDES_LNGTH);  wrChars(rundes);
    date    = rdChars(DATIM_LNGTH);   wrChars(date);
    time    = rdChars(DATIM_LNGTH);   wrChars(time);

    narc = rdI64(); wrI64(narc);
    nsta = rdI64(); wrI64(nsta);

    if nsta ~= 0
        stanam = rdChars(STANAM_LNGTH * double(nsta)); wrChars(stanam);
        idsta  = fread(fin, double(nsta), '*int64');   fwriteIf(fout, idsta, 'int64');
        re     = rdD(1); wrD(re);
        flat   = rdD(1); wrD(flat);
        stacor = rdD(3 * double(nsta)); wrD(stacor);
    end

    % -------------------------
    % 2) Arc Header
    % -------------------------
    ns   = rdI64(); wrI64(ns);
    nf   = rdI64(); wrI64(nf);
    np   = rdI64(); wrI64(np);
    nga  = rdI64(); wrI64(nga);
    iarc = rdI64(); wrI64(iarc);
    utc  = rdD(1);  wrD(utc);

    igrps = rdI64(); wrI64(igrps);

    if igrps > 0
        igrpcnt = fread(fin, double(igrps), '*int64'); fwriteIf(fout, igrpcnt, 'int64');
        igrptot = fread(fin, double(igrps), '*int64'); fwriteIf(fout, igrptot, 'int64');
        %#ok<NASGU>
    end

    idsat  = fread(fin, double(ns), '*int64'); fwriteIf(fout, idsat, 'int64');
    stasat = rdD(6 * double(ns)); wrD(stasat);

    npar = 6*double(ns) + double(np);

    aprval = rdD(npar); wrD(aprval);
    parnam = rdChars(PARNAM_LNGTH * npar); wrChars(parnam);
    parsig = rdD(npar); wrD(parsig);
    parsca = rdD(npar); wrD(parsca);

    % -------------------------
    % 3) Observation records
    % -------------------------
    extra   = double(opts.extraDoubles);
    rowsize = npar + extra + double(igrps);

    assert(all(parCols >= 1 & parCols <= npar), ...
        "parCols must be within 1..npar (npar=%d).", npar);

    if opts.formatCheckOnly
        info = struct('npar',npar,'rowsize',rowsize,'igrps',double(igrps), ...
                      'extra',extra,'note',"Format check only; no records processed.");
        return;
    end

    % Partials start index (1-based) within record vector
    p0 = double(opts.partialsStartInRecord);
    assert(p0 >= 1 && (p0 + npar - 1) <= rowsize, ...
        "partialsStartInRecord makes partials block exceed record length.");

    partialIdxAll = p0 : (p0 + npar - 1);
    overwriteIdx  = partialIdxAll(parCols);   % absolute indices in rec

    % Residual and sigma are in extra block positions 4 and 5 (your correction)
    % => absolute record indices are 4 and 5
    iRes = 4;
    iSig = 5;
    assert(extra >= 5, "extraDoubles must be at least 5 to write residual/sigma at rec(4:5).");

    % newVals must have K+2 columns
    K = numel(parCols);
    assert(size(newVals,2) == K + 2, ...
        "newVals must be Nobs x (K+2): [partials_for_parCols, residual, sigma].");

    nObs = 0;

    while true
        rec = fread(fin, rowsize, '*double');
        if numel(rec) < rowsize
            break; % EOF
        end
        nObs = nObs + 1;

        assert(size(newVals,1) >= nObs, ...
            "newVals has only %d rows but file has at least %d observations.", ...
            size(newVals,1), nObs);

        % 1) overwrite selected partials
        rec(overwriteIdx) = newVals(nObs, 1:K).';

        % 2) overwrite residual and sigma from last two columns
        rec(iRes) = newVals(nObs, K+1);
        rec(iSig) = newVals(nObs, K+2);

        fwrite(fout, rec, 'double');
    end

    info = struct();
    info.npar    = npar;
    info.rowsize = rowsize;
    info.extra   = extra;
    info.igrps   = double(igrps);
    info.nObs    = nObs;
    info.parCols = parCols(:).';
    info.partialsStartInRecord = p0;
    info.residualIndexInRecord = iRes;
    info.sigmaIndexInRecord    = iSig;
end

function tf = fwriteIf(fout, data, precision)
    if fout ~= -1
        fwrite(fout, data, precision);
    end
    tf = true;
end

function cleanupFiles(fin, fout)
    if fin  ~= -1, fclose(fin);  end
    if fout ~= -1, fclose(fout); end
end