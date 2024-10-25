#!/bin/bash

# Bash script to run format_csv.py with appropriate arguments

# Ensure the script is executable by running:
# chmod +x run_format_csv.sh

# Usage:
# ./run_format_csv.sh [input.csv] [output_dir]

# Default input and output values
DEFAULT_INPUT="sorted_deduped_effectors"
DEFAULT_OUTPUT_DIR="../Files/Formatted"

# Check if an input argument is provided, else use default
INPUT_FILE="${1:-$DEFAULT_INPUT}"

# Check if an output directory is provided, else use default
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT_DIR}"

# Run the Python script with the determined input/output paths
python3 format_csv.py --input "$INPUT_FILE" --output_dir "$OUTPUT_DIR"
