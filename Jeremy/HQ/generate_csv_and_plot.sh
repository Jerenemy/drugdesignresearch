#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Save the starting directory
START_DIR=$(pwd)

# Function to run the format CSV script
run_format_csv() {
    echo "Starting CSV formatting..."
    SCRIPT_DIR="$START_DIR/../CSV/Input/Scripts"

    if [ ! -d "$SCRIPT_DIR" ]; then
        echo "Error: Directory '$SCRIPT_DIR' does not exist."
        exit 1
    fi

    cd "$SCRIPT_DIR" || { echo "Failed to cd to $SCRIPT_DIR"; exit 1; }

    if [ -x "./run_format_csv.sh" ]; then
        ./run_format_csv.sh "$@"
        echo "CSV formatting completed."
    else
        echo "Error: run_format_csv.sh is not executable or not found in $SCRIPT_DIR."
        exit 1
    fi

    cd "$START_DIR" > /dev/null
}

# Function to run the evaluate rewards script
run_evaluate_rewards() {
    echo "Starting reward evaluation..."
    SCRIPT_DIR="$START_DIR/../CSV/Output/Scripts"

    if [ ! -d "$SCRIPT_DIR" ]; then
        echo "Error: Directory '$SCRIPT_DIR' does not exist."
        exit 1
    fi

    cd "$SCRIPT_DIR" || { echo "Failed to cd to $SCRIPT_DIR"; exit 1; }

    if [ -x "./run_evaluate_rewards.sh" ]; then
        ./run_evaluate_rewards.sh "$@"
        echo "Reward evaluation completed."
    else
        echo "Error: run_evaluate_rewards.sh is not executable or not found in $SCRIPT_DIR."
        exit 1
    fi

    cd "$START_DIR" > /dev/null
}

# Function to run the generate plot script
run_generate_plot() {
    echo "Starting plot generation..."
    SCRIPT_DIR="$START_DIR/../R/Scripts"

    if [ ! -d "$SCRIPT_DIR" ]; then
        echo "Error: Directory '$SCRIPT_DIR' does not exist."
        exit 1
    fi

    cd "$SCRIPT_DIR" || { echo "Failed to cd to $SCRIPT_DIR"; exit 1; }

    if [ -x "./generate_rewards_plots.sh" ]; then
        ./generate_rewards_plots.sh "$@"
        echo "Plot generation completed."
    else
        echo "Error: generate_rewards_plots.sh is not executable or not found in $SCRIPT_DIR."
        exit 1
    fi

    cd "$START_DIR" > /dev/null
}

# Display usage information
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --format-csv       Run the CSV formatting script."
    echo "  --evaluate-rewards Run the reward evaluation script."
    echo "  --generate-plot    Run the plot generation script."
    echo "  -h, --help         Display this help message."
    echo ""
    echo "If no options are provided, all tasks will be executed in the following order:"
    echo "1. Format CSV"
    echo "2. Evaluate Rewards"
    echo "3. Generate Plot"
}

# Parse command-line arguments
if [ $# -eq 0 ]; then
    # No arguments provided; run all tasks
    run_format_csv
    run_evaluate_rewards
    run_generate_plot
    exit 0
fi

# Iterate over all provided arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --format-csv)
            shift # past argument
            # Pass any subsequent arguments not starting with '--' as parameters to run_format_csv
            FORMAT_CSV_ARGS=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                FORMAT_CSV_ARGS+=("$1")
                shift
            done
            run_format_csv "${FORMAT_CSV_ARGS[@]}"
            ;;
        --evaluate-rewards)
            shift # past argument
            # Pass any subsequent arguments not starting with '--' as parameters to run_evaluate_rewards
            EVAL_REWARDS_ARGS=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                EVAL_REWARDS_ARGS+=("$1")
                shift
            done
            run_evaluate_rewards "${EVAL_REWARDS_ARGS[@]}"
            ;;
        --generate-plot)
            shift # past argument
            # Pass any subsequent arguments not starting with '--' as parameters to run_generate_plot
            GENERATE_PLOT_ARGS=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                GENERATE_PLOT_ARGS+=("$1")
                shift
            done
            run_generate_plot "${GENERATE_PLOT_ARGS[@]}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done
