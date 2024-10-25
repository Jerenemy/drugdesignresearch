#!/bin/bash

# generates the rewards plots with the input file as the csv of the effectors with the 
# rewards calculated, located in ../../CSV/Output/Files/csv_rewards/output_rewards.csv

# call the all_in_1_final.R script:

DEFAULT_CSV="../../CSV/Output/Files/csv_rewards/output_rewards.csv"

# check if an input CSV file provided as arg
if [ "$#" -ge 1 ]; then
    INPUT_CSV="$1"
else 
    INPUT_CSV="$DEFAULT_CSV"
fi

# check if the csv file exists
if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: csv file '$INPUT_CSV' not found!"
    exit 1
fi 

# call the R script with the provided or default input CSV file
Rscript all_in_1_final.R "$INPUT_CSV"

# check if R script executed successfully
if [ $? -eq 0 ]; then   
    echo "plots generated successfully"
else
    echo "error while generating plots"
    exit 1
fi
