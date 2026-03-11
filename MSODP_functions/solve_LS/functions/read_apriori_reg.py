import numpy as np
import struct

def read_reg_apriori(filename: str, endianness: str = "ieee-le"):
    e = endianness.lower()
    ep = "<" if e in ("ieee-le", "little", "le", "<") else ">"

    with open(filename, "rb") as f:
        # --- 1. File Header (Matches C #defines) ---
        f.read(10)  # PROGNAM_LNGTH
        f.read(280) # RUNDES_LNGTH
        f.read(10)  # DATIM_LNGTH (Date)
        f.read(10)  # DATIM_LNGTH (Time)
        
        # Read long long (8 bytes)
        narc = struct.unpack(ep + "q", f.read(8))[0]
        nsta = struct.unpack(ep + "q", f.read(8))[0]
        
        if nsta != 0:
            f.read(10 * nsta) # STANAM_LNGTH * nsta
            f.read(8 * nsta)  # idsta array (long long)
            f.read(8)         # re (double)
            f.read(8)         # flat (double)
            f.read(24 * nsta) # stacor (double * 3 * nsta)

        # --- 2. Arc Header ---
        ns = struct.unpack(ep + "q", f.read(8))[0]
        nf = struct.unpack(ep + "q", f.read(8))[0]
        np_ = struct.unpack(ep + "q", f.read(8))[0]
        nga = struct.unpack(ep + "q", f.read(8))[0]
        iarc = struct.unpack(ep + "q", f.read(8))[0]
        f.read(8) # utc (double)

        igrps = struct.unpack(ep + "q", f.read(8))[0]
        if igrps > 0:
            f.read(8 * igrps) # igrpcnt
            f.read(8 * igrps) # igrptot
            
        f.read(8 * ns)      # idsat (long long)
        f.read(48 * ns)     # stasat (double * 6)

        # --- 3. The Parameter Blocks (Matches C logic) ---
        npar = 6 * ns + np_
        
        # Block 1: APRVAL (Read immediately after Arc Header)
        b = f.read(8 * npar)
        aprval = np.frombuffer(b, dtype=np.dtype(ep + "f8"))
        


    return aprval