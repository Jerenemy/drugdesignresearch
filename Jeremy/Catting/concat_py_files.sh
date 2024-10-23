#!/bin/bash

# Set the starting directory to the current directory or a provided argument
START_DIR=${1:-$(pwd)}
OUTPUT_FILE="output.txt"
MAX_DEPTH=""

# Check if a maximum depth argument is provided
if [ -n "$2" ]; then
    MAX_DEPTH="-maxdepth $2"
fi

# Clear the output file if it already exists
> "$OUTPUT_FILE"

# Function to recursively find and concatenate .py files
concatenate_files() {
    local directory=$1
    # Find .py files up to the specified depth
    find "$directory" $MAX_DEPTH -type f -name "*.py" | while read -r file; do
        # Get the relative directory path
        REL_DIR=$(dirname "$file" | sed "s|^$START_DIR|.|")
        echo "\n%%% DIRECTORY = $REL_DIR %%%" >> "$OUTPUT_FILE"
        echo "$$$ FILENAME = $(basename "$file") $$$" >> "$OUTPUT_FILE"
        cat "$file" >> "$OUTPUT_FILE"
        echo "\n" >> "$OUTPUT_FILE"  # Add an extra newline after each file for readability
    done
}

# Run the function starting from the specified directory
concatenate_files "$START_DIR"

echo "All .py files have been concatenated into $OUTPUT_FILE."
