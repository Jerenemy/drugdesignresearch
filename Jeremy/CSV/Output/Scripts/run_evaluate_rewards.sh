#!/bin/bash

# Bash script to run evaluate_rewards.py with appropriate arguments

# Ensure the script is executable by running:
# chmod +x run_evaluate_rewards.sh

# Usage:
# ./run_evaluate_rewards.sh [input.csv] [output.csv] [REWARDS...]

# Default input and output CSV paths
DEFAULT_INPUT="../../Input/Files/Formatted/sorted_deduped_effectors_formatted.csv"
DEFAULT_OUTPUT="../Files/csv_rewards/output_rewards.csv"

# Default list of all rewards if none are specified
ALL_REWARDS="SIZE SYNTH LIPINSKI QED LogP DOCK"

# Check if any arguments are provided; if not, run with default input/output
if [ "$#" -lt 1 ]; then
    echo "No input file provided. Running with default input and output."
    INPUT_FILE="$DEFAULT_INPUT"
    OUTPUT_FILE="$DEFAULT_OUTPUT"
    REWARDS="$ALL_REWARDS"
    APPEND_FLAG=""
else
    # Get the input file from the first argument
    INPUT_FILE="$1"

    # If a second argument is provided and it's not a reward, use it as the output file
    if [ "$#" -ge 2 ] && [[ ! "$2" =~ ^SIZE|SYNTH|QED ]]; then
        OUTPUT_FILE="$2"
        shift 2  # Shift past the input and output arguments
    else
        OUTPUT_FILE="$DEFAULT_OUTPUT"
        shift 1  # Shift past only the input argument
    fi

    # Any remaining arguments are the rewards
    REWARDS="$@"

    # If no rewards are provided, default to using all rewards
    if [ -z "$REWARDS" ]; then
        REWARDS="$ALL_REWARDS"
    fi

    # Handle default input/output in append mode
    if [ "$OUTPUT_FILE" == "$DEFAULT_OUTPUT" ] && [ -f "$DEFAULT_OUTPUT" ]; then
        # If appending, use the input file as output
        OUTPUT_FILE="$INPUT_FILE"
        APPEND_FLAG="--append"
        echo "Output file already exists. Appending to the same file: $OUTPUT_FILE"
    else
        APPEND_FLAG=""
        echo "Creating new output file: $OUTPUT_FILE"
    fi
fi

# Run the Python script with the determined input/output and rewards
python3 evaluate_rewards.py --input "$INPUT_FILE" --output "$OUTPUT_FILE" $APPEND_FLAG --rewards $REWARDS
