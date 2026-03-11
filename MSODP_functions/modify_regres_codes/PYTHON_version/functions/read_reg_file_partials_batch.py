import numpy as np
import struct
import os

def read_reg_batch(
    filename: str,
    start_obs: int,           # 0-based index of the first observation to read
    batch_size: int,          # Number of observations to read in this batch
    *,
    extra_doubles: int = 5,
    endianness: str = "ieee-le"
):
    """
    Read a specific batch of observations from a binary .reg file.
    """
    # -------- endian ----------
    e = endianness.lower()
    ep = "<" if e in ("ieee-le", "little", "le", "<") else ">"

    # -------- helper for header parsing ----------
    def read_exact(f, n):
        b = f.read(n)
        if len(b) != n: raise EOFError
        return b

    def read_i64(f):
        return struct.unpack(ep + "q", read_exact(f, 8))[0]

    # -------- fixed sizes ----------
    PROGNAM_LNGTH = 10
    RUNDES_LNGTH = 280
    DATIM_LNGTH = 10
    STANAM_LNGTH = 10
    PARNAM_LNGTH = 10

    with open(filename, "rb") as f:
        # ---- 1. Header Navigation (finding rowsize) ----
        f.read(PROGNAM_LNGTH + RUNDES_LNGTH + 2 * DATIM_LNGTH)
        narc = read_i64(f)
        nsta = read_i64(f)

        if nsta != 0:
            f.seek(STANAM_LNGTH * nsta + 8 * nsta + 8 + 8 + 8 * (3 * nsta), 1)

        ns = read_i64(f)
        nf = read_i64(f)
        np_ = read_i64(f)
        f.seek(16 + 8, 1) # nga, iarc + f64_array(1)

        igrps = read_i64(f)
        if igrps > 0:
            f.seek(16 * igrps, 1) # Two i64 arrays

        f.seek(8 * ns + 48 * ns, 1) # idsat + stasat
        npar = 6 * ns + np_

        # Skip parameter arrays and names
        f.seek(8 * npar + (PARNAM_LNGTH * npar) + 16 * npar, 1)

        # ---- 2. Navigate and Read Batch ----
        rowsize = npar + extra_doubles + igrps
        record_bytes = rowsize * 8
        dt = np.dtype(ep + "f8")

        # Jump to the specific start observation
        f.seek(start_obs * record_bytes, 1)

        # Read the batch buffer
        buf = f.read(batch_size * record_bytes)
        if not buf:
            return None
            
        actual_read = len(buf) // record_bytes
        
        # Reshape buffer into (N_batch, Rowsize)
        # Using .copy() ensures the memory is independent and safe for manipulation
        records = np.frombuffer(buf, dtype=dt).reshape(actual_read, rowsize).copy()

        # ---- 3. Extract and Format columns ----
        # Indices: 3=residual, 4=sigma, 5+=partials
        residuals = records[:, 3:4]
        sigmas    = records[:, 4:5]
        partials  = records[:, 5:]

        # Return concatenated result: partials followed by [residual, sigma]
        return np.hstack((partials, residuals, sigmas))