from __future__ import annotations
import struct
import numpy as np
import os
import shutil
from typing import Dict, Optional, Sequence

def reg_update_partials_batch(
    in_file: str,
    start_obs: int,           # 0-based index of the first observation to update
    parCols: Sequence[int],   # 1-based indices (1..npar)
    newVals_batch: np.ndarray, # 2D array: (batch_size, K+2)
    *,
    out_file: Optional[str] = None,
    extra_doubles: int = 5,
    partials_start_in_record: Optional[int] = None,
    endianness: str = 'ieee-le'
) -> Dict[str, object]:
    """
    Python version of reg_update_partials_allRows_v2 optimized for batches.
    """
    # 1. Handle File Setup
    target_file = in_file
    if out_file:
        # If output file doesn't exist, copy the input file as a template
        if not os.path.exists(out_file):
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            shutil.copy2(in_file, out_file)
        target_file = out_file

    # 2. Endianness and Helpers
    e = endianness.lower().strip()
    ep = "<" if e in ("ieee-le", "little", "le", "<") else ">"
    dt_f8 = np.dtype(ep + "f8")

    def read_i64(f): return struct.unpack(ep + "q", f.read(8))[0]

    # 3. Validate Inputs
    parCols = np.asarray(parCols, dtype=int).ravel()
    newVals_batch = np.asarray(newVals_batch, dtype=float)
    if newVals_batch.ndim == 1:
        newVals_batch = newVals_batch.reshape(1, -1)
    
    batch_size = newVals_batch.shape[0]
    K = parCols.size
    p0 = partials_start_in_record if partials_start_in_record is not None else (extra_doubles + 1)

    # 4. Open and Navigate Header
    with open(target_file, "r+b") as f:
        # File Header (fixed lengths)
        f.read(10 + 280 + 10 + 10) # PROGNAM, RUNDES, DATES
        narc = read_i64(f)
        nsta = read_i64(f)
        if nsta != 0:
            f.seek(10 * nsta + 8 * nsta + 8 + 8 + (8 * 3 * nsta), 1)

        # Arc Header
        ns = read_i64(f)
        nf = read_i64(f)
        np_val = read_i64(f)
        f.seek(16 + 8, 1) # nga, iarc, utc
        igrps = read_i64(f)
        if igrps > 0:
            f.seek(16 * igrps, 1) # igrpcnt, igrptot

        f.seek(8 * ns + 48 * ns, 1) # idsat, stasat
        
        npar = int(6 * ns + np_val)
        # Skip aprval, parnam, parsig, parsca
        f.seek(8 * npar + (10 * npar) + 16 * npar, 1)

        # 5. Calculate Record Layout
        rowsize = npar + extra_doubles + igrps
        record_bytes = rowsize * 8
        
        # 6. Navigate to Batch and Update
        f.seek(start_obs * record_bytes, 1)
        
        total_bytes = record_bytes * batch_size
        buf = f.read(total_bytes)
        
        if not buf:
            raise EOFError("Start observation index is beyond end of file.")
            
        actual_batch = len(buf) // record_bytes
        data = np.frombuffer(buf, dtype=dt_f8).copy().reshape(actual_batch, rowsize)

        # Apply Updates (Vectorized)
        # partials start at p0-1 (0-based)
        overwrite_idx = (p0 - 1) + (parCols - 1)
        
        # residual=index 3, sigma=index 4 (MATLAB 4 and 5)
        i_res, i_sig = 3, 4 

        data[:, overwrite_idx] = newVals_batch[:actual_batch, :K]
        data[:, i_res] = newVals_batch[:actual_batch, K]
        data[:, i_sig] = newVals_batch[:actual_batch, K + 1]

        # Write back
        f.seek(-len(buf), 1)
        f.write(data.tobytes(order="C"))

    return {
        "npar": npar,
        "rowsize": rowsize,
        "nObs_updated": actual_batch,
        "start_index": start_obs
    }