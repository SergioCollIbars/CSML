#!/bin/bash

# Define paths
OBSERVATION_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/GG_apriori/"
REGRESS_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/regress_reg/"
OUTPUT_PATH="/Users/sergiocollibars/Desktop/CSML/MSODP_functions/NSM_regress_reg/"
BATCH_SIZE=10

echo "Starting the NSM mapping process..."

# Launch the python script using the variables
python3 main.py "$OBSERVATION_PATH" "$REGRESS_PATH" "$OUTPUT_PATH" "$BATCH_SIZE"

echo "Process completed successfully."
