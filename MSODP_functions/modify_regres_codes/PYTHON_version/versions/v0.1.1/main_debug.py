import re
from functions.modify_reg_NSM import modify_reg_NSM
import numpy as np
from pathlib import Path

"""
MODIFY REGRES FILES WITH NSM

Description: launch the modification of the regres files based on 
            the NSM approach. 
            Read observations for a particular date. 
            Look for regres files in that date. Need all components: 
            XX, XY, XZ, YY, YZ, ZZ. Otherwise skip the date
            Save the modified regres
"""

def main():
    folder_obs = Path("/Users/sergiocollibars/Documents/GG_observations/120by120").resolve()
    folder_reg = Path("/Users/sergiocollibars/Documents/regres_files").resolve()
    folder_out = Path("./folder_out").resolve()
    mask       = np.array([1, 1, 1, 1, 1, 1])
    NS_dir     = np.array([1, 1, 1, 0, 0, 0])
    w_scale    = 1e-3 

    # Ensure output folder exists
    folder_out.mkdir(parents=True, exist_ok=True)

    # Define the regression suffixes we are looking for
    suffixes = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']

    # Iterate through OBS files
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
                # Modify regress files with NSM
                modify_reg_NSM(str(obs_file), day_reg_paths, str(folder_out), mask, NS_dir, w_scale)

if __name__ == "__main__":
    main()