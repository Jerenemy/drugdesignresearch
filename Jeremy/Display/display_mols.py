"""
contains funcs to display mols read from csv file created by create_csv_rewards.py
csv format has Smiles, Nicknames, Reward1, Reward2, ...
csv file created using list of Smiles, list of Nicknames

displays and saves png of >=25 mols WITHOUT names

driver: mol_imgs (commented out)
"""
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

import sys
sys.path.append('..')

print(f"cwd = {os.getcwd()}")

def create_molecule_objects(csv_path: str, smiles_column: str = 'Smiles', method: str = 'pandas', delimiter: str = ';') -> list:
    """
    Read a CSV file, extract SMILES strings, and convert them to RDKit molecule objects.

    Parameters:
    csv_path (str): The path to the CSV file.
    smiles_column (str): The name of the column containing SMILES strings. Default is 'smiles'.
    method (str): The method to read the CSV file. Options are 'pandas' or 'manual'. Default is 'pandas'.
    delimiter (str): The delimiter used in the CSV file. Default is ';'.
    
    Returns:
    list: A list of RDKit molecule objects.
    """
    smiles_strings = []

    if method == 'pandas':
        try:
            # Read the CSV file using pandas with the correct delimiter
            data = pd.read_csv(csv_path, delimiter=delimiter)
            # Display the first few rows to ensure it is read correctly
            print(data.head())
            # Extract the SMILES strings
            smiles_strings = data[smiles_column].dropna().tolist()  # Drop any NaN values
        except Exception as e:
            print(f"An error occurred while reading the file with pandas: {e}")
    elif method == 'manual':
        try:
            # Manually read the CSV file
            with open(csv_path, 'r') as f:
                for line in f.readlines()[1:]:
                    smiles_strings.append(line.split(',')[1])  # Assuming the second column contains SMILES strings
        except Exception as e:
            print(f"An error occurred while reading the file manually: {e}")
    else:
        raise ValueError("Invalid method specified. Use 'pandas' or 'manual'.")

    # Convert SMILES strings to RDKit molecule objects
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_strings if smiles]

    # Check how many molecules were converted successfully
    print(f"Number of molecules converted: {len(molecules)}")

    return molecules

## Example usage with the correct column name and delimiter:
# csv_path = 'TrainingData/qm9.csv'  # Replace with your actual file path
# molecules = create_molecule_objects(csv_path, smiles_column='smiles', method='manual')

# Example usage:
# csv_path = '../Training/TrainingData/SmallDrug.csv'  # Replace with your actual file path
# molecules = create_molecule_objects(csv_path)


# Example usage:
# csv_path = 'chembl_33.csv' 
# molecules = create_molecule_objects(csv_path)

# Function to generate and save molecule images
def generate_molecule_image(
    molecules: list, 
    is_ai_generated: bool = True,
    img_folder: str = None, 
    img_name: str = 'molecules_grid', 
    img_size: tuple = (300, 300), 
    start_mol_disp_ind: int = 0,
    end_mol_disp_ind: int = 25,
    mols_per_row: int = 5
) -> None:
    """
    Generate and save an image grid of molecules in higher quality.

    Parameters:
    molecules (list): A list of RDKit molecule objects.
    is_ai_generated (bool): Flag indicating if the molecules are AI-generated. Default is True.
    img_folder (str): The folder where the image will be saved. If None, it is set based on is_ai_generated.
    img_name (str): The base name for the saved image file. Default is 'molecules_grid'.
    img_size (tuple): The size of each molecule image in the grid. Default is (300, 300).
    start_mol_disp_ind (int): Starting index of displaying molecules. Default 0. Precond: start_mol_disp_ind < end_mol_disp_ind <= len(molecules)
    end_mol_disp_ind (int): Ending ind of displaying molecules. Default 15. Precond: start_mol_disp_ind < end_mol_disp_ind <= len(molecules)
    mols_per_row (int): The number of molecules per row in the grid. Default is 5.
    
    Returns:
    None
    """
    # Set the img_folder based on is_ai_generated if not provided
    if img_folder is None:
        img_folder = 'ai_output' if is_ai_generated else 'training_output'

    # put output all in output folder
    to_img_folder = 'Output/'
    img_folder = to_img_folder + img_folder

    # Generate the image grid with MolsToGridImage
    img = Draw.MolsToGridImage(
        [m for m in molecules if m is not None][start_mol_disp_ind:end_mol_disp_ind], molsPerRow=mols_per_row, subImgSize=img_size, useSVG=False
    )

    # Ensure the output directory exists
    os.makedirs(img_folder, exist_ok=True)

    # Check for the existence of files and determine the new filename
    i = 1
    while True:
        filename = f'{img_name}{i}.png'
        filepath = os.path.join(img_folder, filename)
        if not os.path.exists(filepath):
            break
        i += 1

    # Save the image directly using img.save method
    img.save(filepath)
    print(f"Image saved as {filepath}")

# Example usage:
# molecules = create_molecule_objects('chembl_33.csv')
# generate_molecule_image(molecules)
    
# Generate and save the molecule image
# generate_molecule_image(molecules, img_size=(1000, 1000), mols_per_row=5)