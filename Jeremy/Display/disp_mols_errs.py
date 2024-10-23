"""
functions to convert smiles strings to molecules objects and display mols

most up-to-date display script

accepts csv file created by csv_effectors.py of format:
    Smiles, Nicknames, Chembl #, Binding Affinity, Reward1, Reward2, ...

csv file created using dict of format:
    {<nickname>: [<smiles>, <binding affinity>, <chembl #>]}

when called, displays and saves png of >=25 mols WITH names

supposedly catches errors of invalid smiles strings and only displays valid ones

driver: mol_imgs_aff.py
"""


import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

def create_molecule_objects(csv_path: str, smiles_column: str = 'Smiles', nickname_column: str = 'Nickname', chembl_column: str = 'Chembl #', binding_column: str = 'Binding Affinity', method: str = 'pandas', delimiter: str = ',') -> list:
    """
    Read a CSV file, extract SMILES strings, and convert them to RDKit molecule objects.

    Parameters:
    csv_path (str): The path to the CSV file.
    smiles_column (str): The name of the column containing SMILES strings. Default is 'Smiles'.
    nickname_column (str): The name of the column containing nicknames. Default is 'Nickname'.
    chembl_column (str): The name of the column containing Chembl numbers. Default is 'Chembl #'.
    binding_column (str): The name of the column containing binding affinity. Default is 'Binding Affinity'.
    method (str): The method to read the CSV file. Options are 'pandas' or 'manual'. Default is 'pandas'.
    delimiter (str): The delimiter used in the CSV file. Default is ','.
    
    Returns:
    list: A list of tuples containing RDKit molecule objects and their nicknames.
    """
    smiles_strings = []
    nicknames = []
    chembl_numbers = []
    binding_affinities = []

    if method == 'pandas':
        try:
            # Read the CSV file using pandas with the correct delimiter
            data = pd.read_csv(csv_path, delimiter=delimiter)
            # Display the first few rows to ensure it is read correctly
            # print(data.head())
            # Extract the SMILES strings, nicknames, chembl numbers, and binding affinities
            smiles_strings = data[smiles_column].dropna().tolist()  # Drop any NaN values
            nicknames = data[nickname_column].dropna().tolist()
            chembl_numbers = data[chembl_column].dropna().tolist()
            binding_affinities = data[binding_column].dropna().tolist()
        except Exception as e:
            print(f"An error occurred while reading the file with pandas: {e}")
    elif method == 'manual':
        try:
            # Manually read the CSV file
            with open(csv_path, 'r') as f:
                for line in f.readlines()[1:]:
                    parts = line.split(delimiter)
                    smiles_strings.append(parts[1])  # Assuming the second column contains SMILES strings
                    nicknames.append(parts[2])  # Assuming the third column contains nicknames
                    chembl_numbers.append(parts[3])  # Assuming the fourth column contains Chembl numbers
                    binding_affinities.append(parts[4])  # Assuming the fifth column contains binding affinities
        except Exception as e:
            print(f"An error occurred while reading the file manually: {e}")
    else:
        raise ValueError("Invalid method specified. Use 'pandas' or 'manual'.")


    # Collect only non-None columns and their corresponding data
    columns_data = [
        (smiles_strings, smiles_column),
        (nicknames, nickname_column),
        (chembl_numbers, chembl_column),
        (binding_affinities, binding_column)
    ]
    valid_columns = [(data, name) for data, name in columns_data if data]

    if not valid_columns:
        print("No valid columns to process.")
        return []

    # Ensure all valid columns have the same length
    min_length = min(len(column[0]) for column in valid_columns)
    for i in range(len(valid_columns)):
        valid_columns[i] = (valid_columns[i][0][:min_length], valid_columns[i][1])

    # Convert SMILES strings to RDKit molecule objects, logging any errors
    molecules = []
    for data in zip(*[column[0] for column in valid_columns]):
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append((mol, *data[1:]))
        else:
            print(f"Failed to convert SMILES: {smiles}")
    print(f"Number of molecules converted: {len(molecules)}")

    return molecules

    # # Convert SMILES strings to RDKit molecule objects, logging any errors
    # molecules = []
    # print(smiles_strings)
    # print(zip(smiles_strings, nicknames, chembl_numbers, binding_affinities))
    # for smiles, nickname, chembl, binding in zip(smiles_strings, nicknames, chembl_numbers, binding_affinities):
    #     mol = Chem.MolFromSmiles(smiles)
    #     print(f"mol = {mol}, smiles = {smiles}, nickname= {nickname}, chembl = {chembl}, binding = {binding}")
    #     if mol is not None:
    #         molecules.append((mol, nickname, chembl, binding))
    #     else:
    #         print(f"Failed to convert SMILES: {smiles}")

    # Check how many molecules were converted successfully

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
    molecules (list): A list of tuples containing RDKit molecule objects and their nicknames.
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

    # Put output all in output folder
    to_img_folder = 'Output/'
    img_folder = to_img_folder + img_folder

    # Filter out any molecules that are None, assuming molecules has tuples with at least two values
    mols_to_display = [(m, nickname) for m, nickname, *rest in molecules if m is not None][start_mol_disp_ind:end_mol_disp_ind]

    # Generate the image grid with MolsToGridImage
    valid_mols = []
    valid_legends = []

    for mol, nickname in mols_to_display:
        try:
            img = Draw.MolToImage(mol, size=img_size)
            valid_mols.append(mol)
            valid_legends.append(nickname)
        except Exception as e:
            print(f"An error occurred while processing molecule {nickname}: {e}")

    if valid_mols:
        try:
            img = Draw.MolsToGridImage(
                valid_mols, molsPerRow=mols_per_row, subImgSize=img_size, legends=valid_legends, useSVG=False
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
        except Exception as e:
            print(f"An error occurred while generating the image: {e}")
    else:
        print("No valid molecules to display.")

# Example usage:
# molecules = create_molecule_objects('chembl_33.csv')
# generate_molecule_image(molecules)

# Generate and save the molecule image
# generate_molecule_image(molecules, img_size=(1000, 1000), mols_per_row=5)
