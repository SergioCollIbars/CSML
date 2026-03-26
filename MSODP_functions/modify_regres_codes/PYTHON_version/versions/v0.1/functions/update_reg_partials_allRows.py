from __future__ import annotations
import struct
import numpy as np
from typing import Dict, Optional, Sequence


def reg_update_partials(
    inReg: str,
    outReg: str,
    parCols: Sequence[int],
    newVals: np.ndarray,
    *,
    extra_doubles: int = 5,
    partials_start_in_record: Optional[int] = None, 
    endianness: str = "ieee-le",
    format_check_only: bool = False,
) -> Dict[str, object]:
    """
    Update a binary .reg file by overwriting selected partial columns AND residual+sigma.

    Parameters
    ----------
    inReg, outReg : str
        Input/output .reg files (binary).
    parCols : sequence of int
        1-based indices (1..npar) of the parameter partials to overwrite. (res & sigam are always written)
    newVals : ndarray (Nobs x (K+2))
        Columns:
          - 0..K-1 : partials corresponding to parCols
          - K      : residual  (writes to record element 4 in MATLAB)
          - K+1    : sigma     (writes to record element 5 in MATLAB)
    extra_doubles : int
        Number of doubles at start of each record (must be >= 5 to write residual/sigma).
    partials_start_in_record : int (1-based)
        Where the partials block begins within a record. Default = extra_doubles + 1.
    endianness : {"ieee-le","ieee-be"}
        Little or big endian.
    format_check_only : bool
        If True, only parses header and returns layout info; no output written.

    Returns
    -------
    info : dict
        Metadata: npar, rowsize, nObs, etc.
    """
    # --- endian ---
    e = endianness.lower().strip()
    if e in ("ieee-le", "little", "le", "<"):
        ep = "<"
    elif e in ("ieee-be", "big", "be", ">"):
        ep = ">"
    else:
        raise ValueError("endianness must be 'ieee-le' or 'ieee-be'.")

    # --- validate inputs ---
    parCols = np.asarray(parCols, dtype=int).ravel()
    if parCols.size == 0 or np.any(parCols <= 0):
        raise ValueError("parCols must be non-empty and 1-based positive.")

    newVals = np.asarray(newVals, dtype=float)
    if newVals.ndim != 2:
        raise ValueError("newVals must be 2D (Nobs x (K+2)).")

    if partials_start_in_record is None:
        partials_start_in_record = extra_doubles + 1  # 1-based like MATLAB

    if extra_doubles < 5:
        raise ValueError("extra_doubles must be >= 5 to overwrite residual/sigma at record(4:5).")

    # residual and sigma indexing
    i_res = 3
    i_sig = 4

    # --- small binary helpers ---
    def read_exact(f, n: int) -> bytes:
        b = f.read(n)
        if len(b) != n:
            raise EOFError
        return b

    def read_i64(f) -> int:
        return struct.unpack(ep + "q", read_exact(f, 8))[0]

    def write_i64(f, x: int) -> None:
        f.write(struct.pack(ep + "q", int(x)))

    def read_i64_array(f, n: int) -> np.ndarray:
        b = read_exact(f, 8 * n)
        return np.frombuffer(b, dtype=np.dtype(ep + "i8"), count=n)

    def write_i64_array(f, arr) -> None:
        a = np.asarray(arr, dtype=np.dtype(ep + "i8"))
        f.write(a.tobytes(order="C"))

    def read_f64_array(f, n: int) -> np.ndarray:
        b = read_exact(f, 8 * n)
        return np.frombuffer(b, dtype=np.dtype(ep + "f8"), count=n)

    def write_f64_array(f, arr) -> None:
        a = np.asarray(arr, dtype=np.dtype(ep + "f8"))
        f.write(a.tobytes(order="C"))

    # ---------- fixed header sizes ----------
    PROGNAM_LNGTH = 10
    RUNDES_LNGTH = 280
    DATIM_LNGTH = 10
    STANAM_LNGTH = 10
    PARNAM_LNGTH = 10

    with open(inReg, "rb") as fin:
        fout = None if format_check_only else open(outReg, "wb")

        try:
            # 1) file header
            prognam = read_exact(fin, PROGNAM_LNGTH)
            rundes = read_exact(fin, RUNDES_LNGTH)
            date = read_exact(fin, DATIM_LNGTH)
            time = read_exact(fin, DATIM_LNGTH)

            narc = read_i64(fin)
            nsta = read_i64(fin)

            if fout is not None:
                fout.write(prognam); fout.write(rundes); fout.write(date); fout.write(time)
                write_i64(fout, narc); write_i64(fout, nsta)

            if nsta != 0:
                stanam = read_exact(fin, STANAM_LNGTH * int(nsta))
                idsta = read_i64_array(fin, int(nsta))
                re_ = float(read_f64_array(fin, 1)[0])
                flat = float(read_f64_array(fin, 1)[0])
                stacor = read_f64_array(fin, 3 * int(nsta))

                if fout is not None:
                    fout.write(stanam)
                    write_i64_array(fout, idsta)
                    write_f64_array(fout, [re_])
                    write_f64_array(fout, [flat])
                    write_f64_array(fout, stacor)

            # 2) arc header
            ns = read_i64(fin)
            nf = read_i64(fin)
            np_ = read_i64(fin)
            nga = read_i64(fin)
            iarc = read_i64(fin)
            utc = float(read_f64_array(fin, 1)[0])

            igrps = read_i64(fin)

            if fout is not None:
                write_i64(fout, ns); write_i64(fout, nf); write_i64(fout, np_)
                write_i64(fout, nga); write_i64(fout, iarc)
                write_f64_array(fout, [utc])
                write_i64(fout, igrps)

            if igrps > 0:
                igrpcnt = read_i64_array(fin, int(igrps))
                igrptot = read_i64_array(fin, int(igrps))
                if fout is not None:
                    write_i64_array(fout, igrpcnt)
                    write_i64_array(fout, igrptot)

            idsat = read_i64_array(fin, int(ns))
            stasat = read_f64_array(fin, 6 * int(ns))

            if fout is not None:
                write_i64_array(fout, idsat)
                write_f64_array(fout, stasat)

            npar = 6 * int(ns) + int(np_)

            aprval = read_f64_array(fin, npar)
            parnam = read_exact(fin, PARNAM_LNGTH * npar)
            parsig = read_f64_array(fin, npar)
            parsca = read_f64_array(fin, npar)

            if fout is not None:
                write_f64_array(fout, aprval)
                fout.write(parnam)
                write_f64_array(fout, parsig)
                write_f64_array(fout, parsca)

            # record layout
            rowsize = int(npar + int(extra_doubles) + int(igrps))

            if np.any(parCols > npar):
                raise ValueError(f"parCols must be within 1..npar (npar={npar}).")

            if format_check_only:
                return {
                    "npar": npar,
                    "rowsize": rowsize,
                    "extra_doubles": int(extra_doubles),
                    "igrps": int(igrps),
                    "note": "format_check_only=True, no records processed",
                }

            # partials block indices inside record (0-based)
            p0 = int(partials_start_in_record)  # 1-based
            if not (p0 >= 1 and (p0 + npar - 1) <= rowsize):
                raise ValueError("partials_start_in_record causes partials block to exceed record length.")

            partial_block = np.arange(p0 - 1, (p0 - 1) + npar, dtype=int)
            overwrite_idx = partial_block[parCols - 1]  # parCols are 1-based

            K = int(parCols.size)
            if newVals.shape[1] != K + 2:
                raise ValueError("newVals must have K+2 columns: [partials, residual, sigma].")

            dt_f8 = np.dtype(ep + "f8")
            nObs = 0

            while True:
                buf = fin.read(rowsize * 8)
                if len(buf) == 0:
                    break
                if len(buf) != rowsize * 8:
                    break  # partial record at EOF, stop

                rec = np.frombuffer(buf, dtype=dt_f8, count=rowsize).copy()

                if nObs >= newVals.shape[0]:
                    raise ValueError(
                        f"newVals has {newVals.shape[0]} rows but file has at least {nObs+1} observations."
                    )

                # overwrite selected partials
                rec[overwrite_idx] = newVals[nObs, :K]

                # overwrite residual & sigma (last two cols)
                rec[i_res] = newVals[nObs, K]
                rec[i_sig] = newVals[nObs, K + 1]

                fout.write(rec.tobytes(order="C"))
                nObs += 1

            return {
                "npar": npar,
                "rowsize": rowsize,
                "extra_doubles": int(extra_doubles),
                "igrps": int(igrps),
                "nObs": nObs,
                "parCols": parCols.tolist(),
                "partials_start_in_record": int(p0),
                "residual_index_matlab": 4,
                "sigma_index_matlab": 5,
            }

        finally:
            if fout is not None:
                fout.close()