"""
creates csv_file from LIST of format 
    [(<nickname>, [<smiles>, <binding affinity>, <chembl #>])

UPDATED TO TAKE A LIST OF TUPLES

can also append to that csv file from dict

csv file name: affinity_rewards.csv

called by mol_imgs_aff.py 

has no driver file, run this with 
    python csv_effectors.py

Functions:
- create_csv(effectors_list, reward_module, output_file): 
  Calculates rewards for each effector and writes them to a new CSV file.
- append_to_csv(effectors_list, reward_module, output_file, start_line): 
  Appends new data to an existing CSV file starting from the specified line.
"""

from rdkit import Chem
import csv
# from torch.utils.tensorboard import SummaryWriter
import sys
sys.path.append('..')
sys.path.append('../..')
from Rewards.rewards import FinalRewardModule
from TrainingMain.config_utils import generateRewardModule

def create_csv(effectors_list, reward_module, output_file):
    ## Calculate rewards for each test molecule
    # print(f"hello, {effectors_list}, {output_file}")
    for key, data in effectors_list:
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print(f"\nADDING MOLECULE: {key}...")  # Print nickname
            reward_module.GiveReward(mol)

    # Collect reward data from each SingleReward object
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}
    import os

    # if not os.path.exists(output_file):
    #     os.chdir('CSV')
    try:
        # Write the rewards to a CSV file
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            header = ['Smiles', 'Nickname', 'Chembl #', 'Binding Affinity'] + reward_names
            csvwriter.writerow(header)

            for i, (key, data) in enumerate(effectors_list):
                
                smiles, binding_affinity, chembl_number = data
                row = [smiles, key, chembl_number, binding_affinity] + [reward_data[reward][i] for reward in reward_names]
                csvwriter.writerow(row)
    except Exception as e:
        print(f"Error {e} writing csv, i = {i}, (key, data) = {(key, data)}, len(effectors_list) = {len(effectors_list)}")
    

def append_to_csv(effectors_list, reward_module, output_file, start_line):
    # Calculate rewards for each test molecule
    for key, data in effectors_list:  # UPDATE: Iterate directly over the list of tuples
        smiles = data[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            reward_module.GiveReward(mol)
        else:
            print(f"Failed to convert SMILES: {smiles}")

    # Collect reward data from each SingleReward object
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}

    # Read existing data from the CSV file
    with open(output_file, 'r', newline='') as csvfile:
        existing_data = list(csv.reader(csvfile))

    # Insert new data starting from the specified line
    new_data = []
    for i, (key, data) in enumerate(effectors_list):  # UPDATE: Iterate directly over the list of tuples
        smiles, binding_affinity, chembl_number = data
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            print(f"Appending {key}...")
            row = [smiles, key, chembl_number, binding_affinity] + [reward_data[reward][i] for reward in reward_names]
            new_data.append(row)

    # Combine existing and new data
    updated_data = existing_data[:start_line] + new_data + existing_data[start_line:]

    # Write the combined data back to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for row in updated_data:
            csvwriter.writerow(row)

# Example usage
def main():
    # Create a SummaryWriter instance
    writer = None #SummaryWriter(log_dir='./logs/test')

    # Generate a list of reward modules
    reward_list = generateRewardModule(['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK'], Wandb=False)

    # Initialize the FinalRewardModule
    final_reward_module = FinalRewardModule(writer, reward_list, wandb_log=False, scaling=False)

    # Define the effectors dictionary
    effectors = {
        "Effector1": ["CCn1c2ccccc2c2cc(NC(=O)COC(=O)C3CCC(=O)N3)ccc21", 213000, "CHEMBL1390581"],
        "Effector2": ["CCn1c2ccccc2c2cc(CNC)ccc21", 140000, "CHEMBL1235116"],
        "Effector3": ["O=C(O)c1cc(I)c(-n2cccc2)c(I)c1O", 21000, "CHEMBL4165245"],
        "Effector4": ["CCN(CC)c1nc2c(C(=O)O)c(O)c(I)c(-n3cccc3)c2s1", 4000, "CHEMBL4162051"],
        "Effector5": ["Ic1cc2c([nH][nH]c2=O)c(I)c1-n1cccc1", 50000, "CHEMBL4164831"],
        "Effector6": ["O=C(O)c1cc(-c2ccc(F)cc2)c(-n2cccc2)c(I)c1O", 113000, "CHEMBL4168268"],
        "Effector7": ["CCN(CC)C1CCN(Cc2cc(C#CCNc3ccccc3)cc(I)c2O)CC1", 10000, "CHEMBL4170894"],
        "Effector8": ["COc1ccc(-c2cc(C(=O)O)c(O)c(I)c2-n2cccc2)cc1", 122000, "CHEMBL4171387"],
        "Effector9": ["CCOc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 63000, "CHEMBL4163490"],
        "Effector10": ["CCCOc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 22000, "CHEMBL4168045"],
        "Effector11": ["C=CCOc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 33000, "CHEMBL4160123"],
        "Effector12": ["O=C(O)c1cc(OCCC(F)(F)F)c(-n2cccc2)c(I)c1O", 108000, "CHEMBL4175243"],
        "Effector13": ["O=C(O)c1cc(OCCO)c(-n2cccc2)c(I)c1O", 726000, "CHEMBL4161192"],
        "Effector14": ["O=C(O)c1cc(OCCc2ccccc2)c(-n2cccc2)c(I)c1O", 116000, "CHEMBL4171798"],
        "Effector15": ["CC(C)COc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 30000, "CHEMBL4176302"],
        "Effector16": ["CCCSc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 14000, "CHEMBL4168439"],
        "Effector17": ["CCCCOc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 30000, "CHEMBL4164612"],
        "Effector18": ["O=C(O)c1cc(I)cc(I)c1O", 820000, "CHEMBL1232243"],
        "Effector19": ["O=C(O)c1cc(I)c(-n2cccc2)c(I)c1O", 21000, "CHEMBL4165245"],  # Duplicate name with same value
        "Effector20": ["OCc1cc(I)c(-n2cccc2)c(I)c1O", 20000, "CHEMBL4172014"],
        "Effector21": ["NN=C1CCCc2c1[nH]c1ccc(Cl)cc21", 105000, "CHEMBL2180125"],
        "Effector22": ["O=C(NO)c1cc(I)c(-n2cccc2)c(I)c1O", 139000, "CHEMBL4161418"],
        "Effector23": ["CN(C)Cc1cc(I)c(-n2cccc2)c(I)c1O", 60000, "CHEMBL4175464"],
        "Effector24": ["O=C(O)c1cc(-c2ccccc2)c(-n2cccc2)c(I)c1O", 306000, "CHEMBL4160358"],
        "Effector25": ["CCn1c2ccccc2c2cc(CNC)ccc21", 140000, "CHEMBL1235116"],  # Duplicate name with same value
        "Effector26": ["CN1CCN(C(=O)c2cc(I)c(-n3cccc3)c(I)c2O)CC1", 30000, "CHEMBL4168664"],
        "Effector27": ["CCCCOc1cc(C(=O)O)c(O)c(I)c1-n1cccc1", 30000, "CHEMBL4164612"],  # Duplicate name with same value
        "Effector28": ["CNC(=O)c1cc(I)c(-n2cccc2)c(I)c1O", 70000, "CHEMBL4176523"]
    }

    # Create the CSV file
    create_csv(effectors, final_reward_module, 'files/affinity_rewards.csv')

    # Append new data to the CSV file from the specified line
    new_effectors = {
        "PK11000": ["CSc1ncc(Cl)c(C(=O)O)n1", 26.2, "PK11000"]
    }
    append_to_csv(new_effectors, final_reward_module, 'files/affinity_rewards.csv', start_line=1)

if __name__ == '__main__':
    main()
