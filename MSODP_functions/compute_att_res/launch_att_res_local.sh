#!/bin/bash

# Path with nominal Grav Gradient observations
folder_GG_obs="/Users/sergiocollibars/Documents/GG_observations"

# Path with nominal attitude observations
folder_AttObs_nom="/Users/sergiocollibars/Documents/att_nom"

# Path with true attitude observations
folder_AttObs_true="/Users/sergiocollibars/Documents/att_true"

echo "Starting the Attitude residual anlysis process..."

# Launch the python script using the variables
python3 main.py "$folder_GG_obs" "$folder_AttObs_nom" "$folder_AttObs_true"

echo "Process completed successfully."
