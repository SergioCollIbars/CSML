import re
from pathlib import Path
import argparse
from functions.modify_reg_NSM import modify_reg_NSM

# Define your folders
# folder_obs = Path("/Users/sergiocollibars/Desktop/CSML/MSODP_functions/GG_apriori/")
# folder_reg = Path("/Users/sergiocollibars/Desktop/CSML/MSODP_functions/regress_reg/")
# folder_out = Path("/Users/sergiocollibars/Desktop/CSML/MSODP_functions/NSM_regress_reg")

def main():
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # 2. Add the arguments you want to receive
    # 'path' is the name of the variable, 'help' is what shows up if you type -h
    parser.add_argument("input_obs_folder", help="The path to the observation folder")
    parser.add_argument("input_regress_folder", help="The path to the regress folder")
    parser.add_argument("output_folder", help="The path where results should be saved")
    parser.add_argument("batch_size", help="number of observation read in regress at a time")

    # 3. Parse the arguments from the terminal
    args = parser.parse_args()

    folder_obs = Path(args.input_obs_folder).resolve()
    folder_reg = Path(args.input_regress_folder).resolve()
    folder_out = Path(args.output_folder).resolve()
    batch_size = int(args.batch_size)

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
                modify_reg_NSM(str(obs_file), day_reg_paths, str(folder_out), batch_size)

if __name__ == "__main__":
    main()