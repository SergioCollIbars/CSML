function write_dcm_ck(t_et, C_BI, ckfile, scid, instid, ref)

    % t_et  : 1xN ephemeris time [s past J2000]
    % C_BI  : 3x3xN DCM from inertial frame to body frame
    % ckfile: output CK file, e.g. 'estim_attitude.bc'
    % scid  : spacecraft ID, e.g. -999000
    % instid: CK frame/instrument ID, e.g. -999010
    % ref   : reference frame, usually 'J2000'

    N = length(t_et);

    sclkdp = zeros(1, N);
    quats  = zeros(4, N);
    avvs   = zeros(3, N);

    for k = 1:N

        % ET -> encoded spacecraft clock time
        sclkdp(k) = cspice_sce2c(scid, t_et(k));

        % SPICE quaternion for rotation from ref frame to body frame
        maxIndx = 3 * k; minIndx = maxIndx - 2;
        ROT     = C_BI(minIndx:maxIndx, :);
        quats(:, k) = cspice_m2q(ROT);

        % Angular velocity not included
        avvs(:, k) = [0; 0; 0];

    end

    if exist(ckfile, 'file')
        delete(ckfile);
    end

    handle = cspice_ckopn(ckfile, 'LRO_EST_ATTITUDE_CK', 0);

    begtim = sclkdp(1);

    endtim = sclkdp(end);

    segid  = 'LRO_EST_BUS_ATTITUDE';

    avflag = false;

    % One interpolation interval starting at the first attitude sample.

    % Must be column vector and must coincide with a pointing time.

    starts = sclkdp(1);

    cspice_ckw03( handle, ...
                  begtim, ...
                  endtim, ...
                  instid, ...
                  ref, ...
                  avflag, ...
                  segid, ...
                  sclkdp', ...
                  quats, ...
                  avvs, ...
                  starts );

    cspice_ckcls(handle);

end