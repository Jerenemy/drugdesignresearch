"""
driver file for disp_mols_errs.py

generates images of molecules in csv file created by csv_effectors.py

imports disp_mols_errs (most updated display file version)
calls functions

csv file: affinity_rewards.csv
"""

import sys
sys.path.append('..')
# from display_mols import *


# # Example usage:
# molecules = create_molecule_objects('molecules_rewards.csv', smiles_column='Molecule', delimiter=',')
# print(molecules)
# generate_molecule_image(molecules)

# from Display.display_named_mols import *

# molecules = create_molecule_objects('files/molecules_rewards.csv', delimiter=',')
# generate_molecule_image(molecules)
    
from Display.disp_mols_errs import *

# molecules = create_molecule_objects('files/molecules_rewards.csv', delimiter=',')
# generate_molecule_image(molecules)

def create_imgs(csv_name):
    molecules = create_molecule_objects(f'files/{csv_name}.csv', delimiter=',')
    generate_molecule_image(molecules, img_name=csv_name)


def main(csv_name):
    
    # csv_name = 'affinity_rewards'
    # csv_name = 'bad_rewards'
    # csv_name = 'bad_rewards1'
    # csv_name = 'bad_rewards3'
    # csv_name = 'real_mol_rewards'
    # csv_name = 'control_real_mols'
    
    create_imgs(csv_name)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        input_csv = sys.argv[1]
    else:
        csv_name = 'control_real_mols'

    # output_csv = pick_lines(input_csv, output_csv, file_type, interval, max_lines, start_line)
    # print(output_csv)
    main(input_csv)