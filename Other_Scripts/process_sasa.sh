#!/bin/bash

# Ensure a directory argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash process_pdb_files.sh <directory_with_pdb_files>"
    exit 1
fi

INPUT_DIR="$1"

# Check if directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' not found."
    exit 1
fi

# Define the Python script to use
PYTHON_SCRIPT="sasa.py"

# Loop through all PDB files in the directory
for pdb_file in "$INPUT_DIR"/*.pdb; do
    if [ -f "$pdb_file" ]; then
        echo "Processing $pdb_file ..."
        python3 "$PYTHON_SCRIPT" "$pdb_file" 
    fi
done

echo "Processing complete."
