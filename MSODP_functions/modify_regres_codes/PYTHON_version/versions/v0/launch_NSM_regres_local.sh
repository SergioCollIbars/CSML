#!/bin/bash

# Define paths
OBSERVATION_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/GG_observations"
REGRESS_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/regres_files/"
OUTPUT_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/NSM_regres_files/"

echo "Starting the NSM mapping process..."

# Launch the python script using the variables
python3 main.py "$OBSERVATION_PATH" "$REGRESS_PATH" "$OUTPUT_PATH"

echo "Process completed successfully."
