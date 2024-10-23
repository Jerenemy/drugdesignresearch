from rdkit import Chem
import csv
# from torch.utils.tensorboard import SummaryWriter
import sys
sys.path.append('..')
sys.path.append('../..')
from Rewards.rewards import FinalRewardModule
from TrainingMain.config_utils import generateRewardModule

from Jeremy.CSV.Output.Scripts.evaluate_rewards import create_csv, append_to_csv

def make_csv(csv_name, effectors):
    # Create a SummaryWriter instance
    writer = None #SummaryWriter(log_dir='./logs/test')

    # Generate a list of reward modules
    reward_list = generateRewardModule(['SIZE', 'SYNTH', 'LIPINSKI', 'QED', 'LogP', 'DOCK'], Wandb=False)

    # Initialize the FinalRewardModule
    final_reward_module = FinalRewardModule(writer, reward_list, wandb_log=False, scaling=False)

    # Define the effectors dictionary

    # Create the CSV file
    create_csv(effectors, final_reward_module, f'files/{csv_name}.csv')


import importlib.util

def load_variables(script_path):
    spec = importlib.util.spec_from_file_location("variables", script_path)
    variables = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(variables)
    return variables

def create_tuples_from_variables(variables_module):
    variables = {name: getattr(variables_module, name) for name in dir(variables_module) if not name.startswith("__")}
    values = list(variables.values())
    tuples_list = [(values[i], values[i+1]) for i in range(0, len(values), 2) if i+1 < len(values)]
    return tuples_list


def create_effectors_tup():
    # from sort_effectors import sort_effectors_by_binding_affinity


    script_path = 'csv_tup_lib.py'  # Replace with your script path
    variables_module = load_variables(script_path)
    tuples_list = create_tuples_from_variables(variables_module)
    # for tuple_list in tuples_list:
    #     print(f"\naggg {tuple_list}") 
    return tuples_list

def create_csvs(csv_names):
    name_effectors_tup = create_effectors_tup()
    # print(name_effectors_tup)

    for csv_name in csv_names:

        for effector_name, effectors_tup in name_effectors_tup:
            if csv_name == effector_name:
                print(effector_name)
                print(effectors_tup)

                make_csv(csv_name, effectors_tup)

                print()
                print(f"data <- read.csv('/home/jeremy/RFiles/csvs/{csv_name}.csv')")
                print(f"ggsave(filename = paste(save_dir, '/{csv_name}.png', sep = ''), plot = combined_plot, width = 20, height = 20)")
                

    for csv_name in csv_names:
        print()
        print(f"data <- read.csv('/home/jeremy/RFiles/csvs/{csv_name}.csv')")
        print(f"ggsave(filename = paste(save_dir, '/{csv_name}.png', sep = ''), plot = combined_plot, width = 20, height = 20)")
                

    

def main(csv_names):
    
    # csv_names = ['affinity_rewards', 'bad_rewards3', 'control_real_mols']
    # csv_names = ['more_real_mols']

    create_csvs(csv_names)
    # name_effectors_tup = create_effectors_tup()

    # print(name_effectors_tup)

    



if __name__ == "__main__":
    csv_names = []
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):
            input_csv = sys.argv[i]
            csv_names.append(input_csv)
        print(f"csv_names = {csv_names}")
    else: 
        print('error in csv_driver')
        csv_names = ['more_real_mols']

    
    main(csv_names)