#!/bin/bash

# create_csv.sh

# Save the current directory
current_dir=$(pwd)

# Change to the specified directory
cd /home/jzay/Desktop/GenAITutorial/Training/TrainingData || exit

# Run the Python script in the specified directory and capture the output

# Define the input and output CSV filenames and other variables
input_csv_name=$1
output_csv_name="$input_csv_name"_out
file_type=$2 # "effectors"
# interval=1
# max_lines=26
# max_lines=10
# start_line=100

## Run the Python script in the specified directory and capture the output
## get_lines.py gets lines from a csv_file and   
# output_csv=$(python3 get_lines.py "$input_csv_name" "$output_csv_name" "$file_type" "$interval" "$max_lines" "$start_line") 
# output_csv=$(python3 get_random_lines.py "$input_csv_name" "$output_csv_name" "$file_type" "$max_lines" "$start_line")


# Check if the Python script was successful
if [ $? -eq 0 ]; then
    echo "CSV file created: $output_csv"
else
    echo "Python script failed"
    exit 1
fi

# Copy the created CSV file back to the original directory
cp "${output_csv%.csv}.csv" "$current_dir/Input/CSV_input/Formatted"

# Change back to the original directory
# cd "$current_dir" || exit

# Run more scripts or programs in the current directory
python3 convert_csv_to_list.py $output_csv $file_type
python3 csv_driver.py $output_csv
python3 imgs_driver.py $output_csv

