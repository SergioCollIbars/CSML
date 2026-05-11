import re
from pathlib import Path
import argparse
import numpy as np

from functions.read_reg_file_partials import read_reg


def main():
    parser = argparse.ArgumentParser()

    folder_obs = Path("/Users/sergiocollibars/Documents/GG_observations/120by120")
    folder_reg = Path("/Users/sergiocollibars/Documents/regres_files")

    # Possible regression suffixes
    suffixes = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']

    # Store residuals by component
    residuals_by_comp = {s: [] for s in suffixes}

    # Iterate through OBS files
    for obs_file in sorted(folder_obs.glob("*.ggr")):

        date_match = re.search(r'(\d{4}-\d{2}-\d{2})', obs_file.name)

        if date_match:
            current_date = date_match.group(1)
            print("Processing: " + str(current_date))

            # Look for any available reg files for this date
            for s in suffixes:
                pattern = f"goce_{s}_eggreg_{current_date}_*.reg"
                matches = list(folder_reg.glob(pattern))

                if matches:
                    reg_path = str(matches[0].absolute())

                    # Read regression file
                    H = read_reg(reg_path)

                    # Store residual column
                    residuals_by_comp[s].append(H[:, -1])
                else:
                    print(f"Missing component {s} for date {current_date}")

    # Check if anything was found
    if all(len(v) == 0 for v in residuals_by_comp.values()):
        print("No regression files were found.")
        return

    print("\nRMS values:")

    total_obs = 0

    for s in suffixes:
        if len(residuals_by_comp[s]) > 0:
            # Concatenate all days for this component
            comp_residuals = np.concatenate(residuals_by_comp[s])

            # Compute RMS for this component
            rms = np.sqrt(np.mean(comp_residuals**2))

            total_obs += len(comp_residuals)

            print(f"{s}: {rms:.6e}")
        else:
            print(f"{s}: not available")

    print(f"\nTotal scalar observations processed: {total_obs}")


if __name__ == "__main__":
    main()