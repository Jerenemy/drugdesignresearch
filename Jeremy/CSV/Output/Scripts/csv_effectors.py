#!/usr/bin/env python3
"""
evaluate_rewards.py

This script reads a CSV file containing molecular data, calculates various rewards for each molecule
using the provided reward modules, and writes the results to a new CSV file.

Usage:
    python evaluate_rewards.py --input input.csv --output output_rewards.csv

Dependencies:
    - RDKit
    - Your custom Rewards module
    - Your TrainingMain.config_utils module
    - SA_Score module

Ensure all dependencies are installed and accessible in your Python environment.
"""

import csv
import sys
import argparse
import os
from rdkit import Chem

# Adjust the Python path to include parent directories if necessary
sys.path.append('../../../..')
# sys.path.append('../..')

# Import your custom reward modules
from Rewards.rewards import FinalRewardModule
from TrainingMain.config_utils import generateRewardModule

def create_csv(effectors_list, reward_module, output_file):
    """
    Calculates rewards for each molecule in effectors_list and writes them to output_file.

    Parameters:
        effectors_list (list of tuples): Each tuple contains (Nickname, [Smiles, Binding_Affinity, Chembl#]).
        reward_module (FinalRewardModule): The reward module to calculate rewards.
        output_file (str): Path to the output CSV file.
    """
    # Calculate rewards for each molecule
    for key, data in effectors_list:
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print(f"Adding molecule: {key} ({smiles})")
            reward_module.GiveReward(mol)
        else:
            print(f"Warning: Failed to parse SMILES for {key}: {smiles}")

    # Collect reward data
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}

    # Prepare to write CSV
    header = ['Smiles', 'Nickname', 'Chembl #', 'Binding_Affinity'] + reward_names

    try:
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(header)

            for i, (key, data) in enumerate(effectors_list):
                smiles, binding_affinity, chembl_number = data
                row = [
                    smiles,
                    key,
                    chembl_number,
                    binding_affinity
                ]
                # Append rewards
                for reward in reward_names:
                    try:
                        reward_value = reward_data[reward][i]
                        row.append(reward_value)
                    except IndexError:
                        print(f"Error: Reward index out of range for {reward} on {key}")
                        row.append('')
                csvwriter.writerow(row)
        print(f"Successfully wrote rewards to {output_file}")
    except Exception as e:
        print(f"Error writing to CSV file {output_file}: {e}")

def append_to_csv(effectors_list, reward_module, output_file, start_line):
    """
    Appends new molecular data and their rewards to an existing CSV file.

    Parameters:
        effectors_list (list of tuples): Each tuple contains (Nickname, [Smiles, Binding_Affinity, Chembl#]).
        reward_module (FinalRewardModule): The reward module to calculate rewards.
        output_file (str): Path to the existing CSV file.
        start_line (int): The line number from which to start appending.
    """
    # Calculate rewards for each molecule
    for key, data in effectors_list:
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print(f"Appending molecule: {key} ({smiles})")
            reward_module.GiveReward(mol)
        else:
            print(f"Warning: Failed to parse SMILES for {key}: {smiles}")

    # Collect reward data
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}

    # Read existing data
    try:
        with open(output_file, 'r', newline='') as csvfile:
            existing_data = list(csv.reader(csvfile))
    except FileNotFoundError:
        print(f"Error: The file {output_file} does not exist. Cannot append data.")
        return

    # Prepare new data to append
    new_data = []
    for i, (key, data) in enumerate(effectors_list):
        smiles, binding_affinity, chembl_number = data
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            row = [
                smiles,
                key,
                chembl_number,
                binding_affinity
            ]
            # Append rewards
            for reward in reward_names:
                try:
                    reward_value = reward_data[reward][i]
                    row.append(reward_value)
                except IndexError:
                    print(f"Error: Reward index out of range for {reward} on {key}")
                    row.append('')
            new_data.append(row)

    # Insert new data at the specified start_line
    updated_data = existing_data[:start_line] + new_data + existing_data[start_line:]

    # Write back to the CSV file
    try:
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            for row in updated_data:
                csvwriter.writerow(row)
        print(f"Successfully appended rewards to {output_file}")
    except Exception as e:
        print(f"Error writing to CSV file {output_file}: {e}")

def parse_effectors_dict(effectors_dict):
    """
    Converts a dictionary of effectors to a list of tuples.

    Parameters:
        effectors_dict (dict): Dictionary where key is Nickname and value is [Smiles, Binding_Affinity, Chembl#].

    Returns:
        list of tuples: Each tuple contains (Nickname, [Smiles, Binding_Affinity, Chembl#]).
    """
    return list(effectors_dict.items())

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Calculate molecular rewards and generate a CSV file.")
    parser.add_argument('--input', type=str, required=True, help='Path to the input CSV file containing molecular data.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output CSV file to write rewards.')
    parser.add_argument('--append', action='store_true', help='Flag to append data instead of creating a new CSV.')
    parser.add_argument('--start_line', type=int, default=1, help='Line number to start appending data.')
    args = parser.parse_args()


    input_file = args.input
    output_file = args.output
    append_flag = args.append
    start_line = args.start_line

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: The input file {input_file} does not exist.")
        sys.exit(1)

    # Read the input CSV file
    effectors_dict = {}
    try:
        with open(input_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                nickname = row.get('Nickname')
                smiles = row.get('Smiles')
                binding_affinity = row.get('Binding_Affinity')
                chembl_number = row.get('Chembl#')
                if nickname and smiles:
                    effectors_dict[nickname] = [smiles, binding_affinity, chembl_number]
                else:
                    print(f"Warning: Missing Nickname or Smiles in row: {row}")
    except Exception as e:
        print(f"Error reading input CSV file {input_file}: {e}")
        sys.exit(1)

    effectors_list = parse_effectors_dict(effectors_dict)

    # Generate reward modules
    reward_types = ['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK']
    reward_list = generateRewardModule(reward_types, Wandb=False)

    # Initialize the FinalRewardModule
    writer = None  # Placeholder if you plan to use TensorBoard or similar
    final_reward_module = FinalRewardModule(writer, reward_list, wandb_log=False, scaling=False)

    
    
    # add optionality to specify only selective rewards to calculate
    if append_flag:
        # Append mode
        # add default input/output dir for append_mode (they are now just both the output_dir)
        input_dir = "../Files/csv_rewards/"
        output_dir = "../Files/csv_rewards/"
        append_to_csv(effectors_list, final_reward_module, output_file, start_line=start_line)
    else:
        # Create mode
        # add default input/output dir for create_mode
        input_dir = "../../Input/Files/Formatted/"
        output_dir = "../Files/csv_rewards/"
        create_csv(effectors_list, final_reward_module, output_file)

if __name__ == '__main__':
    main()