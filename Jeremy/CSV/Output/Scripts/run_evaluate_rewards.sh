#!/bin/bash

# Bash script to run evaluate_rewards.py with appropriate arguments

# Ensure the script is executable by running:
# chmod +x run_evaluate_rewards.sh

# Usage:
# ./run_evaluate_rewards.sh [input.csv] [output.csv] [REWARDS...]

# Default file name
DEFAULT_FILE_NAME="sorted_deduped_effectors"
# Default list of all rewards if none are specified
ALL_REWARDS="SIZE SYNTH LIPINSKI QED LogP DOCK"

# Check if any arguments are provided; if not, run with default input/output
if [ "$#" -eq 0 ]; then
    echo "No file name provided. Running with default file name."
    FILE_NAME="$DEFAULT_FILE_NAME"
    REWARDS="$ALL_REWARDS"
    # APPEND_FLAG="--append"
else
    # Get the file name from the first argument

    # If a second argument is provided and it's not a reward, use it as the output file
    if [[ ! $1 =~ ^SIZE|SYNTH|LIPINSKI|QED|LogP|DOCK$ ]]; then
        FILE_NAME="$1"
        # echo "args entered: $@"
        shift  # Shift past the file_name arg
    else
        FILE_NAME="$DEFAULT_FILE_NAME" # no file_name entered, only specific rewards
    fi
    # echo "args entered: $@ not in if"
    # Any remaining arguments are the rewards
    REWARDS="$@"

    # If no rewards are provided, default to using all rewards
    if [ -z "$REWARDS" ]; then
        REWARDS="$ALL_REWARDS"
    fi

    # Handle default input/output in append mode
    if [ -f "../Files/csv_rewards/${FILE_NAME}_rewards.csv" ]; then
        # If appending, use the input file as output
        APPEND_FLAG="--append"
        echo "Output file already exists. Appending to the same file: $OUTPUT_FILE"
    else
        APPEND_FLAG=""
        echo "Creating new output file: $OUTPUT_FILE"
    fi
fi

# Run the Python script with the determined input/output and rewards
# python3 evaluate_rewards.py --input "$INPUT_FILE" --output "$OUTPUT_FILE" $APPEND_FLAG --rewards $REWARDS
python3 evaluate_rewards.py --file_name "$FILE_NAME" $APPEND_FLAG --rewards $REWARDS
