import numpy as np
from pathlib import Path
from functions.convert_quaternion_2_DCM import quat2dcm
from datetime import datetime, timezone

def read_att_to_DCM(filename, t0, tf):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"Could not find file: {filename}")

    print(f"    Loading: {filename}")

    # --- J2000 epoch (UTC) ---
    tJ2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)

    num_records = 0

    with open(path, 'r') as f:
        # Parse header
        for line in f:
            if "+datarecor" in line:
                num_records = int(line.split()[-1])

            if "+eoh______" in line or "+eoh" in line:
                break

        # Read only numeric records
        data = np.loadtxt(f, max_rows=num_records)

    if data.size == 0:
        return np.array([])

    data = np.atleast_2d(data)

    # Columns from file
    years   = data[:, 0].astype(int)
    doys    = data[:, 1].astype(int)
    seconds = data[:, 2].astype(float)

    # Compute seconds past J2000 row by row
    tJ2K = np.empty(len(data), dtype=float)

    for i in range(len(data)):
        base_date = datetime.strptime(f"{years[i]}-{doys[i]}", "%Y-%j").replace(tzinfo=timezone.utc)
        tJ2K[i] = (base_date - tJ2000).total_seconds() + seconds[i]

    # Keep only rows in [t0, tf]
    mask = (tJ2K >= t0) & (tJ2K <= tf)

    if not np.any(mask):
        return np.array([])

    data = data[mask]

    # Extract quaternions
    q1 = data[:, 4]
    q2 = data[:, 5]
    q3 = data[:, 6]
    q0 = data[:, 7]
    q_all = np.column_stack((q0, q1, q2, q3))

    # Convert to DCM
    n = q_all.shape[0]
    DCM_stack = np.zeros((n, 3, 3))
    for k in range(n):
        DCM_stack[k] = quat2dcm(q_all[k])

    return DCM_stack