#!/bin/bash
# DESCRIPTION: read all files with .reg in the current folder and use dmpreg routine to create the regress file
    # Create output directory if it does not exist
    mkdir -p regress_txt

    # Loop over all .reg files in current directory
    for file in *.reg; do

        # Skip if no .reg files are found
        [ -e "$file" ] || continue

        # Extract filename without extension
        base="${file%.reg}"

        # Run dmpreg and redirect output
        dmpreg "$file" 3 -1 > "${base}.txt"

        echo "Processed: $file -> ${base}.txt"
    done 

