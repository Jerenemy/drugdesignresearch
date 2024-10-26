#!/bin/bash

# generates the rewards plots with the input file as the csv of the effectors with the 
# rewards calculated, located in ../../CSV/Output/Files/csv_rewards/output_rewards.csv

# call the all_in_1_final.R script:
DEFAULT_DIR="../../CSV/Output/Files/csv_rewards/"
DEFAULT_CSV_NAME="sorted_deduped_effectors"

# check if an input CSV file provided as arg
if [ "$#" -ge 1 ]; then
    CSV__NAME="$1"
else 
    CSV_NAME="$DEFAULT_CSV_NAME"
fi

CSV_FILE_PATH="${DEFAULT_DIR}${CSV_NAME}_rewards.csv"

# check if the csv file exists
if [ ! -f "$CSV_FILE_PATH" ]; then
    echo "Error: csv file '$CSV_FILE_PATH' not found!"
    exit 1
fi 

# call the R script with the provided or default input CSV file
Rscript all_in_1_final.R "$CSV_FILE_PATH"

# check if R script executed successfully
if [ $? -eq 0 ]; then   
    echo "plots generated successfully"
else
    echo "error while generating plots"
    exit 1
fi
