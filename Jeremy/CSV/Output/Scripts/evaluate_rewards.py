#!/usr/bin/env python3
"""
evaluate_rewards.py

This script reads a CSV file containing molecular data, calculates various rewards for each molecule
using the provided reward modules, and writes the results to a new CSV file or appends to an existing one.

Usage:
    # To create a new CSV with all rewards
    python evaluate_rewards.py --input path/to/input.csv --output path/to/output_rewards.csv

    # To create a new CSV with selective rewards
    python evaluate_rewards.py --input path/to/input.csv --output path/to/output_rewards.csv --rewards SIZE SYNTH QED

    # To append data to an existing CSV with all rewards
    python evaluate_rewards.py --append --output path/to/output_rewards.csv

    # To append data with selective rewards
    python evaluate_rewards.py --append --output path/to/output_rewards.csv --rewards SIZE QED

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
sys.path.append('../../../..')  # Adjust this path based on your project structure
# sys.path.append('../..')  # Uncomment or modify if needed

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
            # print(f"Adding molecule: {key} ({smiles})")
            reward_module.GiveReward(mol)
        else:
            print(f"Warning: Failed to parse SMILES for {key}: {smiles}")

    # Collect reward data
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}

    # Prepare to write CSV
    header = ['Smiles', 'Nickname', 'Chembl#', 'Binding_Affinity'] + reward_names

    try:
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(header)
            print(f"# effectors = {len(effectors_list)-1}")
            for i, (key, data) in enumerate(effectors_list):
                smiles, binding_affinity, chembl_number = data
                row = [
                    smiles,
                    key,
                    chembl_number,
                    binding_affinity
                ]
                # Append rewards
                # print(f"len = {len(reward_names)}")
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
    Updates rewards for existing molecular data in the CSV file. If a molecule does not exist, appends it.

    Parameters:
        effectors_list (list of tuples): Each tuple contains (Nickname, [Smiles, Binding_Affinity, Chembl#]).
        reward_module (FinalRewardModule): The reward module to calculate rewards.
        output_file (str): Path to the existing CSV file.
        start_line (int): The line number from which to start updating/appending data.
    """
    # Calculate rewards for each molecule
    for key, data in effectors_list:
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # print(f"Updating molecule: {key} ({smiles})")
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
            header = existing_data[0]  # Preserve header
            existing_data = existing_data[1:]  # Skip header
    except FileNotFoundError:
        print(f"Error: The file {output_file} does not exist. Cannot update data.")
        return
    except Exception as e:
        print(f"Error reading the existing CSV file {output_file}: {e}")
        return

    # Create a dictionary to index existing molecules by Nickname or SMILES
    existing_dict = {}  # UPDATE: Dictionary to track existing molecules by Nickname or SMILES
    for row in existing_data:
        smiles = row[0]
        nickname = row[1]
        existing_dict[nickname] = row

    # Prepare new/updated data
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

        # UPDATE: Replace row if it exists, otherwise add new
        if key in existing_dict:
            existing_dict[key] = row  # Update existing entry
        else:
            existing_dict[key] = row  # Add new entry

    # Write updated data back to the CSV
    try:
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            # Write the header
            csvwriter.writerow(header)
            # Write updated rows
            for row in existing_dict.values():
                csvwriter.writerow(row)
        print(f"Successfully updated rewards in {output_file}")
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
    parser.add_argument('--input', type=str, default=None, help='Path to the input CSV file containing molecular data.')
    parser.add_argument('--output', type=str, default=None, help='Path to the output CSV file to write rewards.')
    parser.add_argument('--append', action='store_true', help='Flag to append data instead of creating a new CSV.')
    parser.add_argument('--start_line', type=int, default=1, help='Line number to start appending data (0-based index).')
    parser.add_argument('--rewards', type=str, nargs='+', default=None, help='List of rewards to calculate (e.g., SIZE SYNTH QED). If not specified, all rewards are calculated.')
    args = parser.parse_args()

    append_flag = args.append

    # Set default directories based on mode
    if append_flag:
        default_input = "../Files/csv_rewards/output_rewards.csv"
        default_output = "../Files/csv_rewards/output_rewards.csv"
    else:
        default_input = "../../Input/Files/Formatted/sorted_deduped_effectors_formatted.csv"
        default_output = "../Files/csv_rewards/output_rewards.csv"

    # Assign input and output files
    input_file = args.input if args.input else default_input
    output_file = args.output if args.output else default_output

    # Check if input file exists
    if not os.path.exists(input_file):
        # print()
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
    if args.rewards:
        reward_types = args.rewards
    else:
        reward_types = ['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK']
    
    print(f"Selected Rewards: {', '.join(reward_types)}")

    reward_list = generateRewardModule(reward_types, Wandb=False)

    # Initialize the FinalRewardModule
    writer = None  # Placeholder if you plan to use TensorBoard or similar
    final_reward_module = FinalRewardModule(writer, reward_list, wandb_log=False, scaling=False)

    if append_flag:
        # Append mode
        append_to_csv(effectors_list, final_reward_module, output_file, start_line=args.start_line)
    else:
        # Create mode
        create_csv(effectors_list, final_reward_module, output_file)

if __name__ == '__main__':
    main()
