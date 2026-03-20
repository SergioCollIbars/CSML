import re
from pathlib import Path
import numpy as np
from uniplot import plot
import argparse

from functions.GG_rotation_partials import GG_rotation_partials
from functions.compute_att_error import compute_att_error
from functions.read_GG_obs import read_GG_obs

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("folder_GG_obs", help="The path to the gravity gradient observations data folder (no noise)")
    parser.add_argument("folder_AttObs_nom", help="The path to the S/C atittude data folder (no noise)")
    parser.add_argument("folder_AttObs_true", help="The path to the s/C attitude data folder (noise)")

    args = parser.parse_args()

    folder_GG_obs      = Path(args.folder_GG_obs)
    folder_AttObs_nom  = Path(args.folder_AttObs_nom)
    folder_AttObs_true = Path(args.folder_AttObs_true)

    blocks = []
    blocks_2 = []
    blocks_3 = []
    # Iterate through OBS files
    for obs_file in sorted(folder_GG_obs.glob("*.ggr")):
        # Extract date from filename (e.g. 2008-08-01)
        date_match = re.search(r"(\d{4}-\d{2}-\d{2})", obs_file.name)
        if not date_match:
            print(f"Skipping {obs_file.name}: no date found in filename.")
            continue

        current_date = date_match.group(1)
        print(f"Processing: {current_date}")

        pattern = re.compile(rf"grace\.{re.escape(current_date)}_A_61\.att\.pre$")

        # Look for the nominal attitude observation file
        matches_AttObs_nom = None
        for filename in sorted(folder_AttObs_nom.glob("*.pre")):
            if pattern.fullmatch(filename.name):
                matches_AttObs_nom = str(filename)
                break

        # Look for the true attitude observation file
        matches_AttObs_true = None
        for filename in sorted(folder_AttObs_true.glob("*.pre")):
            if pattern.fullmatch(filename.name):
                matches_AttObs_true = str(filename)
                break

        # Skip if either file is missing
        if matches_AttObs_nom is None or matches_AttObs_true is None:
            print(f"Skipping {current_date}: missing attitude file(s).")
            print(f"  nominal: {matches_AttObs_nom}")
            print(f"  true   : {matches_AttObs_true}")
            continue

        # Read Grav. gradient observations for the current date
        G = read_GG_obs(str(obs_file))

        if G.size == 0:
            print(f"Skipping {current_date}: empty GG observation file.")
            continue

        t0 = G[0, 0]
        tf = G[-1, 0]

        # Compute attitude error for the current date
        att_err = compute_att_error(matches_AttObs_nom, matches_AttObs_true, t0, tf)

        # Check consistency between GG observations and attitude errors
        N_att = att_err.shape[1]
        N_gg = G.shape[0]
        N = min(N_att, N_gg)

        if N_att != N_gg:
            print(f"Warning for {current_date}:")
            print(f"  att_err has {N_att} samples")
            print(f"  G has {N_gg} samples")
            print(f"  Using N = {N}")

        # Compute and store GG residual due to 1st-order attitude errors
        dY = np.zeros((N, 6))
        tr = np.zeros((N, 1))
        at = np.zeros((N, 3))

        for k in range(N):
            yaw_err = att_err[2, k]
            pitch_err = att_err[1, k]
            roll_err = att_err[0, k]

            at[k, :] = np.array([yaw_err, pitch_err, roll_err])

            obs = G[k, 1:7]   # 6 GG components
            H = GG_rotation_partials(obs)   # expected shape (6, 3)

            dY[k, :] = H @ np.array([yaw_err, pitch_err, roll_err])
            tr[k, 0] = np.abs(dY[k, 0] + dY[k, 3] + dY[k, 5])

        # Append as (N x 6) block
        blocks.append(dY)

        # Append trace
        blocks_2.append(tr)

        # Append attitude errors
        blocks_3.append(at)

    # Stack all days into one big array
    if blocks:
        big_array = np.vstack(blocks)
        trc_array = np.vstack(blocks_2)
        att_array = np.vstack(blocks_3)
    else:
        big_array = np.empty((0, 6))
        trc_array = np.empty((0, 1))
        att_array = np.empty((0, 3))

    print("Final shape of big_array:", big_array.shape)

    # Compute RMS for each of the 6 columns
    if big_array.size > 0:
        rms_values = np.sqrt(np.mean(big_array**2, axis=0))

        print("RMS values for the 6 GG components:")
        for i, val in enumerate(rms_values, start=1):
            print(f"Component {i}: {val:.6e}")
    else:
        print("big_array is empty, RMS cannot be computed.")

    # save file
    np.savetxt(
        "GG_attitude_residuals.txt",
        big_array,
        fmt="%.12e", 
        delimiter=" "
    )
    np.savetxt(
    "attitude_residuals.txt",
    att_array,
    fmt="%.12e", 
    delimiter=" "
    )
    
    print(f"maximum trace value = {np.max(trc_array)}")

    # plot residuals
    for k in range(0, 6):
        plot(big_array[:, k], y_as_log=False, title=f"Residuals {k+1}")

# run main
if __name__ == "__main__":
    main()