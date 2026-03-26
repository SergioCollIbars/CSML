#!/bin/bash

# Define paths
OBSERVATION_PATH="/Users/sergiocollibars/Documents/GG_observations/120by120"
REGRESS_PATH="/Users/sergiocollibars/Documents/regres_files"
OUTPUT_PATH="./folder_out"

# MASK: XX XY XZ YY YZ ZZ. EX: 111111 (all) / 100101 (diag only)
# NS DIRECTIONS: 1-6.      EX:345 / 4 (0-index based)

echo "Starting the NSM mapping process..."

# Launch the python script using the variables
python3 main.py "$OBSERVATION_PATH" "$REGRESS_PATH" "$OUTPUT_PATH" 111110 00011

echo "Process completed successfully."
