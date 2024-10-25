#!/bin/bash

# contains functionality to 
# 1. generate the formatted csv
# 2. generate/append to the csv with the evaluated rewards
# 3. generate the plot with R

# first set it up to pass no additional args to each of the sh scripts, just cd into their directories and call them
# cd into Input to create formatted csv
cd ../CSV/Input/Scripts || exit 1

# example call: ./generate_csv_and_plot.sh <format csv?> <args that format_csv takes (up to 2)> <evaluate rewards?> <args that eval rewards takes (up to 8 including reward options)> <generate plot?> <args that generate plot takes (up to 1)>

if [ $1 -eq "run_format_csv" ]; then
    if [ $2 -ne "run_evaluate_rewards" || $2 -ne "generate_rewards_plots" ] # and $3 also not them
    ./run_format_csv.sh
    # otherwise call it with the correct cmd line args



cd ../../Output/Scripts || exit 1

# same thing for here
if [ $1 -eq "run_format_csv" ]; then
    if [ $2 -ne "run_evaluate_rewards" || $2 -ne "generate_rewards_plots" ] # and $3 also not them
    ./run_format_csv.sh
    # otherwise call it with the correct cmd line args


cd ../../../R/Scripts || exit 1

# same thing for here