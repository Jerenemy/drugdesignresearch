import csv
import sys
import argparse
import os

def main(input_file, output_dir):

    output_file = f'{output_dir}/{os.path.basename(input_file)}_formatted.csv'
    input_file = f"../Files/Raw/{input_file}.csv"

    # Open the input and output files
    with open(input_file, 'r', newline='', encoding='utf-8') as csv_infile, \
        open(output_file, 'w', newline='', encoding='utf-8') as csv_outfile:

        reader = csv.DictReader(csv_infile)
        fieldnames = ['Smiles', 'Nickname', 'Chembl#', 'Binding_Affinity', 'Affinity_Units', 
                      'SizeReward', 'SynthReward', 'LipinskiReward', 'QED_Reward', 'LogPReward', 'DockReward']
        writer = csv.DictWriter(csv_outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            output_row = {
                'Smiles': row['Ligand SMILES'],
                'Nickname': row['Unnamed: 0'],
                'Chembl#': '',
                'Binding_Affinity': '',
                'Affinity_Units': '',
                'SizeReward': '',
                'SynthReward': '',
                'LipinskiReward': '',
                'QED_Reward': '',
                'LogPReward': '',
                'DockReward': ''
            }

            # Find the first non-empty affinity value
            if row['Ki (nM)']:
                output_row['Binding_Affinity'] = row['Ki (nM)']
                output_row['Affinity_Units'] = 'Ki (nM)'
            elif row['IC50 (nM)']:
                output_row['Binding_Affinity'] = row['IC50 (nM)']
                output_row['Affinity_Units'] = 'IC50 (nM)'
            elif row['Kd (nM)']:
                output_row['Binding_Affinity'] = row['Kd (nM)']
                output_row['Affinity_Units'] = 'Kd (nM)'
            elif row['EC50 (nM)']:
                output_row['Binding_Affinity'] = row['EC50 (nM)']
                output_row['Affinity_Units'] = 'EC50 (nM)'

            writer.writerow(output_row)
        

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Format CSV file with molecular data.")
    parser.add_argument('--input', type=str, default='sorted_deduped_effectors', help='Input CSV file name without extension.')
    parser.add_argument('--output_dir', type=str, default='../Files/Formatted', help='Directory to save the formatted CSV.')

    args = parser.parse_args()

    main(args.input, args.output_dir)
