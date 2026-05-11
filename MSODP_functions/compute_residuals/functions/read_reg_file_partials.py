import numpy as np
import struct
import os


def read_reg(
    filename: str,
    *,
    extra_doubles: int = 5,
    endianness: str = "ieee-le"
):
    """
    Read a binary .reg file and return ALL observation record columns.

    Returns
    -------
    records  -> (Nobs, rowsize)
    """
    reg_name = os.path.basename(filename)
    print(f"Loading: {reg_name}")
    
    # -------- endian ----------
    e = endianness.lower()
    if e in ("ieee-le", "little", "le", "<"):
        ep = "<"
    elif e in ("ieee-be", "big", "be", ">"):
        ep = ">"
    else:
        raise ValueError("endianness must be 'ieee-le' or 'ieee-be'.")

    # -------- helpers ----------
    def read_exact(f, n):
        b = f.read(n)
        if len(b) != n:
            raise EOFError
        return b

    def read_i64(f):
        return struct.unpack(ep + "q", read_exact(f, 8))[0]

    def read_i64_array(f, n):
        b = read_exact(f, 8 * n)
        return np.frombuffer(b, dtype=np.dtype(ep + "i8"))

    def read_f64_array(f, n):
        b = read_exact(f, 8 * n)
        return np.frombuffer(b, dtype=np.dtype(ep + "f8"))

    # -------- fixed sizes ----------
    PROGNAM_LNGTH = 10
    RUNDES_LNGTH = 280
    DATIM_LNGTH = 10
    STANAM_LNGTH = 10
    PARNAM_LNGTH = 10

    with open(filename, "rb") as f:

        # ---- file header ----
        f.read(PROGNAM_LNGTH)
        f.read(RUNDES_LNGTH)
        f.read(DATIM_LNGTH)
        f.read(DATIM_LNGTH)

        narc = read_i64(f)
        nsta = read_i64(f)

        if nsta != 0:
            f.read(STANAM_LNGTH * nsta)
            read_i64_array(f, nsta)
            read_f64_array(f, 1)
            read_f64_array(f, 1)
            read_f64_array(f, 3 * nsta)

        # ---- arc header ----
        ns = read_i64(f)
        nf = read_i64(f)
        np_ = read_i64(f)
        nga = read_i64(f)
        iarc = read_i64(f)
        read_f64_array(f, 1)

        igrps = read_i64(f)

        if igrps > 0:
            read_i64_array(f, igrps)
            read_i64_array(f, igrps)

        read_i64_array(f, ns)
        read_f64_array(f, 6 * ns)

        npar = 6 * ns + np_

        read_f64_array(f, npar)
        f.read(PARNAM_LNGTH * npar)
        read_f64_array(f, npar)
        read_f64_array(f, npar)

        # ---- record layout ----
        rowsize = npar + extra_doubles + igrps
        dt = np.dtype(ep + "f8")

        records = []

        while True:
            buf = f.read(rowsize * 8)
            if len(buf) == 0:
                break
            if len(buf) != rowsize * 8:
                break

            rec = np.frombuffer(buf, dtype=dt)
            residual = rec[3]
            records.append(residual)

        records = np.vstack(records)

    return records