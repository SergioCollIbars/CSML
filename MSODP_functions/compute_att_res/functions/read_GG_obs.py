import re
import numpy as np
from datetime import datetime, timezone
import os


def read_GG_obs(file_path: str) -> np.ndarray:
    """
    READ Grav. Grad (GG) observations

    Parameters
    ----------
    file_path : str
        Path to the .txt/.ggr file. The filename contains the date 'YYYY-MM-DD'.

    Returns
    -------
    out : ndarray, shape (N, 7)
        Columns: [time_day_sec, data(:,21:26)]
        where time_day_sec = data(:,3) + secJ2K and secJ2K is seconds past J2000.
    """
    if not file_path:
        print("Warning: No file provided.")
        return np.empty((0, 7), dtype=float)

    # --- regex date in filename ---
    pat = r"(?P<date>\d{4}-\d{2}-\d{2})"
    m = re.search(pat, file_path)
    if m is None:
        raise ValueError(f"Could not find YYYY-MM-DD date in filename: {file_path}")

    dstr = m.group("date")

    # --- J2000 epoch (UTC) ---
    t0 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    t = datetime.strptime(dstr, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    secJ2K = (t - t0).total_seconds()
    
    file_name = os.path.basename(file_path)
    print(f"    Loading: {file_name}")

    # --- read first line to infer number of columns ---
    with open(file_path, "r") as f:
        first_line = f.readline()
        if first_line == "":
            return np.empty((0, 7), dtype=float)
        nCols = len(first_line.strip().split())

    # --- load all floats (like textscan %f) ---
    raw = np.loadtxt(file_path, dtype=float)  # may be 2D already if well-formed

    if raw.ndim == 1:
        if raw.size % nCols != 0:
            raise ValueError(
                f"File has {raw.size} numbers, not divisible by inferred nCols={nCols}."
            )
        data = raw.reshape((-1, nCols))
    else:
        data = raw
        if data.shape[1] != nCols:
            
            flat = data.ravel()
            if flat.size % nCols != 0:
                raise ValueError(
                    f"Data size {flat.size} not divisible by inferred nCols={nCols}."
                )
            data = flat.reshape((-1, nCols))

    # time_day_sec = data[]:,2] + secJ2K
    time_day_sec = data[:, 2] + secJ2K

    # data(:,20:25) 
    gg_cols = data[:, 20:26]  # 6 columns

    out = np.column_stack([time_day_sec, gg_cols])  # (N, 7)

    # find(~isnan(dataOutput(:,1)));
    mask = ~np.isnan(out[:, 0])
    out = out[mask, :]

    return out