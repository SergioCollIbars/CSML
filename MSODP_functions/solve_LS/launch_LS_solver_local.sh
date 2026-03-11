#!/bin/bash

# Define paths
REGRES_PATH="/Users/sergiocollibars/Documents/regres_2008"

echo "Starting the LS solver  process..."

# Launch the python script using the variables
python3 main.py "$REGRES_PATH"

echo "Process completed successfully."
