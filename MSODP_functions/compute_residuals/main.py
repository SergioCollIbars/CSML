import re
from pathlib import Path
import argparse
import numpy as np

from functions.read_reg_file_partials import read_reg

def main():
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # 2. Add the arguments you want to receive
    # 'path' is the name of the variable, 'help' is what shows up if you type -h
    parser.add_argument("input_obs_folder", help="The path to the observation folder")
    parser.add_argument("input_regress_folder", help="The path to the regress folder")

    # 3. Parse the arguments from the terminal
    args = parser.parse_args()

    folder_obs = Path(args.input_obs_folder).resolve()
    folder_reg = Path(args.input_regress_folder).resolve()

    # Define the regression suffixes we are looking for
    suffixes = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']

    # Iterate through OBS files
    blocks = []
    for obs_file in sorted(folder_obs.glob("*.ggr")):  # Adjust extension if needed
        
        # Extract date from filename (e.g., 2008-08-01)
        # This regex looks for 4 digits-2 digits-2 digits
        date_match = re.search(r'(\d{4}-\d{2}-\d{2})', obs_file.name)
        if date_match:
            current_date = date_match.group(1)
            day_reg_paths = []
            print("Processing: " + str(current_date))
            
            # Look for the 6 specific reg files for this date
            found_all = True
            for s in suffixes:
                # Pattern matches: goce_XX_eggreg_2008-08-01...
                pattern = f"goce_{s}_eggreg_{current_date}_*.reg"
                matches = list(folder_reg.glob(pattern))
                
                if matches:
                    # Convert Path object to a string and append to list
                    day_reg_paths.append(str(matches[0].absolute()))
                else:
                    print(f"Missing component {s} for date {current_date}")
                    found_all = False
                    break
            
            # Process only if we have the full set of 6
            if found_all:
                # compute residuals
                H_XX = read_reg(day_reg_paths[0])
                H_XY = read_reg(day_reg_paths[1])
                H_XZ = read_reg(day_reg_paths[2])
                H_YY = read_reg(day_reg_paths[3])
                H_YZ = read_reg(day_reg_paths[4])
                H_ZZ = read_reg(day_reg_paths[5])

                r = np.array([H_XX[:,-2],H_XY[:, -2],H_XZ[:, -2],\
                    H_YY[:, -2],H_YZ[:, -2],H_ZZ[:, -2] ])
                
                blocks.append(r)
                
    if len(blocks) == 0:
        print("No complete set of regression files was found.")
        return

    # stack all days into one big array of shape (6, Ntotal)
    big_array = np.hstack(blocks)
    print(f"Number of obs. processed: {len(big_array.T)}")

    # compute RMS value for each of the 6 observation channels
    rms = np.sqrt(np.mean(big_array**2, axis=1))

    # print results
    print("\nRMS values:")
    for comp, val in zip(suffixes, rms):
        print(f"{comp}: {val:.6e}")

if __name__ == "__main__":
    main()