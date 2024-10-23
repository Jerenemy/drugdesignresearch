"""
has no driver file, run this with 
    python create_csv_rewards.py

creates a csv file
csv format has Smiles, Nicknames, Reward1, Reward2, ...
csv file created using list of Smiles, list of Nicknames

needs to be called by display_named_mols.py
csv file name: molecules_rewards.csv

different format from affinity_rewards.csv
(does not contain Binding Affinity, Chembl# columns) 

Example csv file format:

Smiles,Nickname,SizeReward,SynthReward,LipinskiReward,QED Reward,LogPReward,DockReward
CSc1ncc(Cl)c(C(=O)O)n1,pk11000,-1.5,-2.0,2.0,1.355334081124468,-1.0,-2
COc1ccc(-c2cc(C(=O)O)c(O)c(I)c2-n2cccc2)cc1,effector 8,1.5,-1.1799614859252676,2.0,1.4128404917045687,1.0,-2

"""

from rdkit import Chem
import sys
import csv
# from torch.utils.tensorboard import SummaryWriter
sys.path.append('..')
sys.path.append('../..')
from Rewards.rewards import FinalRewardModule, SynthReward, QEDReward
from TrainingMain.config_utils import generateRewardModule

def create_csv(smiles_list, nicknames_list, reward_module, output_file):
    
    # Calculate rewards for each test molecule
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        reward_module.GiveReward(mol)

    # Collect reward data from each SingleReward object
    reward_names = [reward.name() for reward in reward_module.r_list]
    reward_data = {reward.name(): reward.reward_list for reward in reward_module.r_list}

    # Write the rewards to a CSV file
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        header = ['Smiles', 'Nickname'] + reward_names
        csvwriter.writerow(header)

        for i, smiles in enumerate(smiles_list):
            row = [smiles, nicknames_list[i]] + [reward_data[reward][i] for reward in reward_names]
            csvwriter.writerow(row)

# Example usage
def main():
    # Create a SummaryWriter instance
    writer = None #SummaryWriter(log_dir='./logs/test')

    # Generate a list of reward modules
    reward_list = generateRewardModule(['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK'], Wandb=False)

    # Initialize the FinalRewardModule
    final_reward_module = FinalRewardModule(writer, reward_list, wandb_log=False, scaling=False)

    # Create a list of test molecules
    smiles_list = [
        'CSc1ncc(Cl)c(C(=O)O)n1', # pk11000
        'COc1ccc(-c2cc(C(=O)O)c(O)c(I)c2-n2cccc2)cc1', # effector 8: binding affinity = -9.6
        'OCc1cc(I)c(-n2cccc2)c(I)c1O', # effector 20: binding affinity (kcal/mol) = -7.4
        'O=C(NO)c1cc(I)c(-n2cccc2)c(I)c1O', # effector 22: binding affinity = -8.9
        "O=C(N[C@H]1COC1=O)c1ccccc1", # random mol downloaded from training file
        'Cc1c2c3c(c(NPCl)c1PC2)CN=N3', # theos test mols v
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'C1CCCCC1',
        'c1ccccc1',
        'CCO'
    ]

    # Create a list of nicknames for the test molecules
    nicknames_list = [
        'pk11000',
        'effector 8',
        'effector 20',
        'effector 22',
        'random real mol',
        'Molecule6',
        'Molecule7',
        'Molecule8',
        'Molecule9',
        'Molecule10'
    ]

    # Create the CSV file
    create_csv(smiles_list, nicknames_list, final_reward_module, 'files/molecules_rewards.csv')

if __name__ == '__main__':
    main()
