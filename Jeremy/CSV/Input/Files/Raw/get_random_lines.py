"""
This script RANDOMLY selects a specified number of lines from an input CSV file and writes them to 
an output CSV file, starting from a given line number. It handles different file types ('qm9' or 
'chembl') by adjusting the delimiter and 'SMILES' column index accordingly. The script can be run 
with command-line arguments or default values.

Must be called in directory /home/jzay/Desktop/GenAITutorial/Training/TrainingData

Called in create_csv.sh: 
    output_csv=$(python3 get_random_lines.py "$input_csv_name" "$output_csv_name" "$file_type" "$max_lines" "$start_line")

Example usage:
    python3 get_random_lines.py SmallDrug new_chembl_csv chembl 26 1

Prints (outputs) the name of the csv_file it created, so it can be captured in the bash script.
"""
import csv
import random
from typing import List

def pick_random_lines(input_csv: str, output_csv: str, file_type: str = 'qm9', max_lines: int = 26, start_line: int = 0) -> None:
    """
    Randomly selects lines from an input CSV file and writes them to an output CSV file.

    Parameters:
        input_csv (str): The path to the input CSV file.
        output_csv (str): The path to the output CSV file.
        file_type (str): The type of file ('qm9' or 'chembl'). Determines the delimiter and the index of the 'smiles' column.
        max_lines (int): The maximum number of lines to write to the output CSV.
        start_line (int): The line number to start processing from in the input CSV.

    Returns:
        None
    """
    delimiter = ',' if file_type == 'qm9' else ';'
    
    with open(input_csv, mode='r', encoding='utf-8-sig') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = list(reader)

    header = lines[0]
    data_lines = lines[start_line+1:]  # Exclude the header and start line
    
    selected_lines: List[List[str]] = [header]
    
    if len(data_lines) < max_lines:
        selected_indices = range(len(data_lines))
    else:
        selected_indices = random.sample(range(len(data_lines)), max_lines)
    
    for index in selected_indices:
        row = data_lines[index]
        smiles_index = header.index('smiles') if file_type == 'qm9' else header.index('Smiles')
        smiles = row[smiles_index]
        if smiles:  # Check if smiles string is present
            selected_lines.append(row)
    
    with open(output_csv, mode='w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter=delimiter)
        writer.writerows(selected_lines)

def main(argv: List[str]) -> None:
    """
    Main function to execute the pick_random_lines function with command-line arguments or default values.

    Parameters:
        argv (List[str]): List of command-line arguments.

    Returns:
        None
    """
    if len(argv) > 1:
        input_csv_name = argv[1]
        output_csv_name = argv[2]
        file_type = argv[3]
        max_lines = int(argv[4])
        start_line = int(argv[5])
    else:
        input_csv_name = 'chembl_1'  # Replace with your input CSV file path
        output_csv_name = 'chembl_25'  # Replace with your desired output CSV file path
        file_type = 'chembl'  # Specify 'qm9' or 'chembl'
        max_lines = 26
        start_line = 0  # Specify the line number to start processing from

    input_csv_path = input_csv_name + '.csv'
    output_csv_path = output_csv_name + '.csv'
    
    # Create the CSV file
    pick_random_lines(input_csv_path, output_csv_path, file_type, max_lines=max_lines, start_line=start_line)

    # Output the name of the created CSV file
    print(output_csv_name)

if __name__ == "__main__":
    import sys
    main(sys.argv)
