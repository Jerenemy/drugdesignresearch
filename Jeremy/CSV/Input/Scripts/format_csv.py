import csv
import sys
def main(input_file, output_dir):

    # input_file = 'input.csv'   # Replace with your input file name
    output_file = f'{output_dir}/{input_file}_out.csv' # Replace with your desired output file name
    input_dir = "../Files/Raw"
    input_file = f"{input_dir}/{input_file}.csv"

    # Open the input and output files
    with open(input_file, 'r', newline='', encoding='utf-8') as csv_infile, \
        open(output_file, 'w', newline='', encoding='utf-8') as csv_outfile:

        reader = csv.DictReader(csv_infile)
        # Define the output fieldnames
        fieldnames = ['Smiles', 'Nickname', 'Chembl#', 'Binding_Affinity', 'Affinity_Units', 
                    'SizeReward', 'SynthReward', 'LipinskiReward', 'QED_Reward', 'LogPReward', 'DockReward']
        writer = csv.DictWriter(csv_outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            output_row = {}
            output_row['Smiles'] = row['Ligand SMILES']
            # output_row['Nickname'] = row['BindingDB Reactant_set_id']
            output_row['Nickname'] = row['Unnamed: 0']
            output_row['Chembl#'] = ''  # Leave blank or populate as needed
            # Find the first non-empty affinity value
            affinity_value = ''
            affinity_units = ''
            if row['Ki (nM)']:
                affinity_value = row['Ki (nM)']
                affinity_units = 'Ki (nM)'
            elif row['IC50 (nM)']:
                affinity_value = row['IC50 (nM)']
                affinity_units = 'IC50 (nM)'
            elif row['Kd (nM)']:
                affinity_value = row['Kd (nM)']
                affinity_units = 'Kd (nM)'
            elif row['EC50 (nM)']:
                affinity_value = row['EC50 (nM)']
                affinity_units = 'EC50 (nM)'

            output_row['Binding_Affinity'] = affinity_value
            output_row['Affinity_Units'] = affinity_units
            # Leave the reward columns blank
            output_row['SizeReward'] = ''
            output_row['SynthReward'] = ''
            output_row['LipinskiReward'] = ''
            output_row['QED_Reward'] = ''
            output_row['LogPReward'] = ''
            output_row['DockReward'] = ''

            writer.writerow(output_row)


if __name__ == "__main__":
    # csv_names = []
    if len(sys.argv) > 1:
        # for i in range(1, len(sys.argv)):
        input_file = sys.argv[1]
    else: 
        input_file = 'sorted_deduplicated_data_effectors'
        
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else: 
        output_dir = '../Files/Formatted'
            # print(i)
            # csv_names.append(input_csv)
        # print(f"csv_names = {csv_names}")
        # print('error in csv_driver')
        # csv_names = ['more_real_mols']

    
    main(input_file, output_dir)