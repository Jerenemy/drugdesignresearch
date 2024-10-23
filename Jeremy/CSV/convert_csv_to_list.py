"""
This script reads an input CSV file and converts its data into a specific format suitable for another script. 
It processes the CSV file based on the specified file type ('qm9' or 'chembl') and generates a tuple list, 
which is then written to a Python file as a variable. The script can be run with command-line arguments or 
default values.

outputs 2 vars into .py file.

Example output:
    var_name = 'var_name'
    var_name_tup = [
        ('Molecule 1', ['CCN', -1, 'CHEMBL1']),
    ]

var_name is a cmd-line arg and gives       

tup format:
    (nickname, [smiles, binding affinity, molecule #])

"""


# convert_csv_to_list.py
import csv
import sys

def generate_effectors_tup(input_csv, file_type='qm9'):
    effectors_tup = []
    delimiter = ',' if file_type == 'qm9' else ';'
    
    with open(input_csv, mode='r', encoding='utf-8-sig') as infile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        print(f"\n\n\nfile type = {file_type}, type = {type(file_type)}\n\n\n")
        for row in reader:
            if file_type == 'qm9':
                mol_id = row['mol_id']
                smiles = row['smiles']
                name = mol_id  # Using mol_id as name for qm9 format
            elif file_type == 'chembl':
                mol_id = row['Parent Molecule']
                smiles = row['Smiles']
                name = row['Name']
            effectors_tup.append((name, [smiles, -1, mol_id]))
    
    return effectors_tup

def write_to_py_file(var_name, effectors_tup, output_py_file):
    with open(output_py_file, mode='a') as pyfile:
        pyfile.write("\n")
        pyfile.write(f"{var_name} = '{var_name}'\n")


        pyfile.write("\n")
        pyfile.write(f"{var_name}_tup = [\n")
        for effector in effectors_tup:
            pyfile.write(f"    {effector},\n")
        pyfile.write("]\n")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        csv_name = sys.argv[1]
    else:
        csv_name = "chembl_25"
    if len(sys.argv) > 2:
        file_type = sys.argv[2]
    else:
        file_type = 'chembl'  # Specify 'qm9' or 'chembl']
        
    input_csv = f'files/{csv_name}.csv'  # Replace with your input CSV file path
    output_py_file = 'csv_tup_lib.py'  # Specify the output .py file path
    var_name = csv_name  # Specify the variable name

    write_to_py_file(var_name, generate_effectors_tup(input_csv, file_type=file_type), output_py_file)
