from __future__ import annotations
import struct
import numpy as np
from typing import Dict, Optional, Sequence


def reg_update_partials(
    reg_file: str,
    obs_index: int,           # 0-based index of the observation to update
    parCols: Sequence[int],
    newVals: np.ndarray,      # 1D array: [partial1, partial2, ..., residual, sigma]
    *,
    extra_doubles: int = 5,
    partials_start_in_record: Optional[int] = None,
    endianness: str = "ieee-le",
) -> Dict[str, object]:
    """
    Updates a SINGLE observation in a .reg file in-place.
    """
    e = endianness.lower().strip()
    ep = "<" if e in ("ieee-le", "little", "le", "<") else ">"

    # --- validate inputs ---
    parCols = np.asarray(parCols, dtype=int).ravel()
    newVals = np.asarray(newVals, dtype=float).ravel() # Expecting 1D for single obs
    
    K = parCols.size
    if newVals.size != K + 2:
        raise ValueError(f"newVals must have {K+2} elements (partials + residual + sigma).")

    if partials_start_in_record is None:
        partials_start_in_record = extra_doubles + 1

    # Open in r+b for in-place modification
    with open(reg_file, "r+b") as f:
        # Helper to skip bytes
        def skip(n): f.seek(n, 1)
        def read_i64(): return struct.unpack(ep + "q", f.read(8))[0]
        def read_f64(): return struct.unpack(ep + "f", f.read(8))[0]

        # 1) Skip File Header
        f.read(10 + 280 + 10 + 10) # PROGNAM, RUNDES, DATES
        narc = read_i64()
        nsta = read_i64()

        if nsta != 0:
            skip(10 * nsta + 8 * nsta + 8 + 8 + (8 * 3 * nsta))

        # 2) Skip Arc Header
        ns = read_i64()
        nf = read_i64()
        np_ = read_i64()
        nga = read_i64()
        iarc = read_i64()
        skip(8) # utc
        igrps = read_i64()

        if igrps > 0:
            skip(8 * igrps * 2)

        skip(8 * ns + 8 * 6 * ns) # idsat, stasat
        npar = 6 * int(ns) + int(np_)

        # Skip aprval, parnam, parsig, parsca
        skip(8 * npar + (10 * npar) + 8 * npar + 8 * npar)

        # 3) Calculate Record Position
        rowsize = npar + extra_doubles + igrps
        record_bytes = rowsize * 8
        
        # Move to the specific observation
        f.seek(obs_index * record_bytes, 1)
        
        # Read the current record to modify it
        buf = f.read(record_bytes)
        if not buf or len(buf) != record_bytes:
            raise IndexError(f"Observation index {obs_index} is out of bounds.")
        
        rec = np.frombuffer(buf, dtype=np.dtype(ep + "f8")).copy()

        # 4) Overwrite values in the array
        p0 = int(partials_start_in_record)
        partial_block_start = p0 - 1
        overwrite_idx = partial_block_start + (parCols - 1)

        rec[overwrite_idx] = newVals[:K]
        rec[3] = newVals[K]     # residual
        rec[4] = newVals[K + 1] # sigma

        # 5) Write back to the SAME location
        f.seek(-record_bytes, 1) # Backtrack to start of this record
        f.write(rec.tobytes(order="C"))

    return {"status": "success", "obs_updated": obs_index, "npar": npar}