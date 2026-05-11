#!/bin/bash

# Define paths & weight
REGRESS_PATH="/Users/sergiocollibars/Documents/regres_files"
WEIGHT=1E-5

echo "Starting the weigth writting process..."

# Launch the python script using the variables
python3 main.py "$REGRESS_PATH" "$WEIGHT"

echo "Process completed successfully."
