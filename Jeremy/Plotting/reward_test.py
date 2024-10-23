"""
file to test rewards from ..Rewards.rewards

not in use

creates a plot showing graph of rewards for each mol using plt
now i create the csv and export it and create the plot with RStudio

"""

import torch
from torch.utils.tensorboard import SummaryWriter
from rdkit import Chem
import sys
sys.path.append('..')
from Rewards.rewards import FinalRewardModule, SynthReward, QEDReward
from TrainingMain.config_utils import generateRewardModule

# def generateRewardModule(reward_names, Wandb=True):
#     rewards = []
#     if "SYNTH" in reward_names:
#         rewards.append(SynthReward(Wandb=Wandb))
#     if "QED" in reward_names:
#         rewards.append(QEDReward(Wandb=Wandb))
#     return rewards

def main():
    # Create a SummaryWriter instance
    writer = SummaryWriter(log_dir='./logs/test')

    # Generate a list of reward modules
    reward_list = generateRewardModule(['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK'], Wandb=False)

    # reward_list = generateRewardModule(['SYNTH', 'QED'], Wandb=False)

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

    # Calculate rewards for each test molecule
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)

        reward = final_reward_module.GiveReward(mol)
        print(f"Reward for molecule {smiles}: {reward}")

    # Print the rewards
    final_reward_module.print_rewards()
    print("plotting")
    # Plot the rewards
    final_reward_module.PlotRewards()

    print("plotted")
if __name__ == '__main__':
        main()