from abc import ABC, abstractmethod
import heapq
import numpy as np
import sys
sys.path.append('../../../..')  # Adjust this path based on your project structure

if __name__ == '__main__':
    from autoDock import VinaWrapper
else:
    from .autoDock import VinaWrapper

# import rdkit
from rdkit import Chem
from rdkit.Chem import Lipinski, Draw
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import qed
import matplotlib.pyplot as plt
import wandb


# import torch
# import math


if __name__ == '__main__':
    from SA_Score import sascorer
else:
    from .SA_Score import sascorer



# lss wandb

class SingleReward(ABC):
    @abstractmethod
    def __init__(self, Wandb = True):
        self.reward_list = []
        self.Wandb = Wandb
        if Wandb:
            wandb.define_metric(self.name(), step_metric="num_mols")
    @abstractmethod
    def _reward(self,mol):
        pass
    
    @abstractmethod
    def name(self):
        pass
    
    def rescale(self, x):
        return x
    
    def giveReward(self, x):
        reward = self.rescale(self._reward(x))
        if self.Wandb:
            wandb.log({self.name(): reward})
        self.reward_list.append(reward)
        return reward
    
    @staticmethod
    def gaussianReward(x, ideal, sigma):
        """
        Calculate the reward based on Gaussian distribution.
        
        Parameters:
        x (float): The value to evaluate.
        ideal (float): The ideal value.
        sigma (float): The standard deviation of the Gaussian distribution.
        
        Returns:
        float: The reward in range (0, 1]
        """
        return np.exp(-0.5 * ((x - ideal) / sigma) ** 2)

        # # Example usage
        # ideal_value = 100
        # sigma = 10  # This determines the "width" of the ideal range
        # value_to_evaluate = 105

        # reward = gaussian_reward(value_to_evaluate, ideal_value, sigma)

    @staticmethod
    def scale_to_range(value, orig_range=(0,1), new_range=(-2,2)):
        """
        Scale a single value from an original range to a new range.
        
        Parameters:
        value (float): The value to be scaled.
        orig_range (tuple): The original range (min, max) of the value.
        new_range (tuple): The new range (min, max) for the value.
        
        Returns:
        float: The scaled value.
        """
        orig_min, orig_max = orig_range
        new_min, new_max = new_range
        
        try:
            # Scale the value to the new range
            scaled_value = new_min + (value - orig_min) * (new_max - new_min) / (orig_max - orig_min)
        except ZeroDivisionError:
            print(f"Error: Division by zero encountered while scaling value {value}. Returning original value.")
            return value
        
        return scaled_value

    @staticmethod
    def transform(score, threshold=0.5, emphasize_above=True, steepness=10):
        """
        Apply a sigmoid transformation to the score to emphasize values
        above or below a certain threshold.
        
        Parameters:
        score (float): The score to be transformed. Expected range is (0,1).
        threshold (float): The threshold value to emphasize around. Expected range is (0,1). (default is 0.5)
        emphasize_above (bool): Whether to emphasize values above the threshold (default is True).
        steepness (float): The steepness of the sigmoid function (default is 10). Expected range is any positive real number.
        
        Returns:
        float: The transformed reward, in the range [0, 1].
        
        Example:
        >>> Scaler.transformed_reward(0.3, threshold=0.5, emphasize_above=True)
        0.18242552380635635
        >>> Scaler.transformed_reward(0.7, threshold=0.5, emphasize_above=False)
        0.18242552380635635
        """
        try:
            if emphasize_above:
                return 1 / (1 + np.exp(-steepness * (score - threshold)))
            else:
                return 1 / (1 + np.exp(steepness * (score - threshold)))
        except ZeroDivisionError:
            print(f"Error: Division by zero encountered while transforming reward for score {score}. Returning score as reward.")
            return score
        
    @staticmethod
    def normalize_value(value, orig_range):
        """
        Normalize a value to the range [0, 1] based on the original range.
        
        Parameters:
        value (float): The value to be normalized.
        orig_range (tuple): The original range (min, max) of the value.
        
        Returns:
        float: The normalized value in the range [0, 1].
        """
        orig_min, orig_max = orig_range
        
        try:
            normalized_value = (value - orig_min) / (orig_max - orig_min)
        except ZeroDivisionError:
            print(f"Error: Division by zero encountered while normalizing value {value}. Returning 0.")
            return 0
        
        return normalized_value
    
    @staticmethod
    def combined_reward(x, transition_point, ideal, sigma):
        """
        Calculate the reward based on a combined linear and Gaussian distribution.
        
        Parameters:
        x (float): The value to evaluate.
        transition_point (float): The point at which the distribution transitions from linear to Gaussian.
        ideal (float): The ideal value.
        sigma (float): The standard deviation of the Gaussian distribution.
        
        Returns:
        float: The reward.
        """
        if x <= transition_point:
            # Linear distribution up to the transition point
            return (x / transition_point) * np.exp(-0.5 * ((x - ideal) / sigma) ** 2)
        else:
            # Gaussian distribution after the transition point
            return np.exp(-0.5 * ((x - ideal) / sigma) ** 2)


class SynthReward(SingleReward):
    def __init__(self, Wandb = True):
        super(SynthReward,self).__init__(Wandb)
        try:
            sascorer.calculateScore(Chem.MolFromSmiles("CC"))
            print('Synth Reward Succesfully initialized')
        except:
            print('synth fail')
            raise 

    def name(self):
        return 'SynthReward'
    
    
    def rescale(self, r):
        """
        returns range [-2,2]
        """
        # divide r by 9, since sascorer.calculateScore gives in range [1,9] (1 is best)
        norm_r = SingleReward.normalize_value(r, (1,10))
        # transform r to emphasize vals above 0.5 (after normalized)
        transf_r = SingleReward.transform(norm_r, emphasize_above=False, threshold=0.37, steepness=10)
        return SingleReward.scale_to_range(transf_r, orig_range=(0,1), new_range=(-2,2))



    
    def _reward(self,mol: Chem.rdchem.Mol):
        """
        returns synthetic accessibility score in range [1,10] (1 is best)
        """
        # print(Chem.MolToSmiles(mol))
        try:
            sascore = sascorer.calculateScore(mol)
            # print(f"sascore = {sascore}")
        except ZeroDivisionError:
            sascore = 10
        return sascore
    
        
    
    
import multiprocessing
# import signal
# import time

# Other imports and class definitions...

class DockReward(SingleReward):
    def __init__(self, receptor_path, Wandb=True):
        super(DockReward, self).__init__(Wandb)
        
        # for Jeremy: if calling rewards.py from different location, relative receptor_path
        #  will be different. therefore check that path exists, if it doesn't, assume it's
        #  being called from Jeremy/CSV/Output/Scripts instead of the usual place, and update 
        #  relative receptor_path accordingly.
        import os
        if not os.path.exists(receptor_path):
            print(f"RECEPTOR_PATH = {receptor_path}")
            receptor_path = '../../../' + receptor_path
        if os.path.exists(receptor_path):
            print("running from Jeremy/CSV/Output/Scripts. receptor path exists")
        else: 
            print(f"invalid receptor_path. using hard-coded receptor path. cwd = {os.getcwd()}")
            receptor_path = '/home/jzay/Desktop/GoodThesis1/DrugDesignThesis/CLEAN/Rewards/y220c_av.pdbqt'
            # raise Exception
        
        self.vinaWrapper = VinaWrapper(receptor_path)
        self.timeout_occurred = None

    def name(self):
        # print("dock")
        return "DockReward"

    def rescale(self, value):
        # print("dock", value)
        # if value > 2:
        #     if self.timeout_occurred:
        #         print("timeout")
        #         self.timeout_occurred = False
        #         return -1.5
        #     print(f"value={value}")
        #     return -2
        reward = (-value) - 5
        return np.clip(reward, a_min=-2, a_max=10)

        # if value > 2:
        #     if self.timeout_occurred:
        #         self.timeout_occurred = False
        #         return -1.5 * 2
        #     return -2 * 2
        # reward = ((-value) - 5)*2
        # return np.clip(reward, a_min=-4, a_max=20)

    def _reward(self, mol, timeout_seconds=20):
        smile = Chem.MolToSmiles(mol)
        energies = self.vinaWrapper.CalculateEnergies(smile)
        print(f"smile = {smile}, energies = {energies}")
        # print(f"receptor_path = {self.receptor_path}")
        import os
        # print(f"os.listdir = {os.listdir('../../../../')}")
        return energies
        # timeout_seconds = 1

        # def calculate_energies(queue, smile):
        #     energies = self.vinaWrapper.CalculateEnergies(smile)
        #     queue.put(energies)

        # queue = multiprocessing.Queue()
        # process = multiprocessing.Process(target=calculate_energies, args=(queue, smile))
        # process.start()
        # process.join(timeout_seconds)

        # if process.is_alive():
        #     # print(f"Timeout: Docking took longer than {timeout_seconds} seconds, terminating process.")
        #     process.terminate()
        #     process.join()
        #     self.timeout_occurred = True
        #     return 3  # or any high penalty value indicating a timeout

        # energies = queue.get() if not queue.empty() else 3
        # print(f"energies = {energies}")
        # return energies

class SizeReward(SingleReward):
    def __init__(self, Wandb = True):
        super(SizeReward,self).__init__(Wandb)
        
    def name(self):
        return "SizeReward"
    
    def rescale(self,value):
        """
        chosen 20.52 because average of num mols in all of marks molecules (the "best" molecules we have)
        4 = std of best mols they came up with, increasing to allow more variability
        
        translate to range (-2,2] from (0,1]
        """
        # target_atoms = 20.52
        # distribution = 10
        # turning_point = 15
        # gauss_r = SingleReward.combined_reward(value, turning_point, target_atoms, distribution)
        # # gauss_r = SingleReward.gaussianReward(value, target_atoms, distribution)
        # return SingleReward.scale_to_range(gauss_r)
        
        target_atoms = 24
        distribution = 10
        transition_point = 20
        gauss_r = SingleReward.combined_reward(value, transition_point, target_atoms, distribution)
        return SingleReward.scale_to_range(gauss_r) #,new_range=(-10,10))
        
    def _reward(self,mol):
        """
        Returns num atoms in molecule.
        """
        # print(f"mol.GetNumAtoms() = {mol.GetNumAtoms()}")
        return mol.GetNumAtoms() # for now, more atoms is "better"

        
class LogPReward(SingleReward):
    """
    Reward for how hydrophopic/hydrophilic molecule is (lipophilic/lipophobic).
    
    LogP > 0 ==> more hydrophobic (lipophilic): less soluble, breaks through cell membranes better.
    LogP < 0 ==> more hydrophilic (lipophobic): more soluble, doesn't break through cell membranes as well.

    Tuned to prefer slightly hydrophobic (+3.3).
    """
    def __init__(self, Wandb = True):
        super(LogPReward,self).__init__(Wandb)
        
    def name(self):
        return 'LogPReward'
    
    def rescale(self, x):
        """
        3.32138 = avg(LogP) of best mols
        0.7 = std of best mols

        translate to range (-2,2] from (0,1]
        """
        gauss_r = SingleReward.gaussianReward(x, 3.94, 1)
        return SingleReward.scale_to_range(gauss_r)


    def _reward(self, mol):
        # print(f"LogP = {MolLogP(mol)}")
        return MolLogP(mol)

    
class LipinskiReward(SingleReward):
    """
    Lipinski's Rule of 5 for Hydrogen Bind Donors (HBD) and Hydrogen Bond Acceptors (HBA).

    Rule: compound should have no more than 5 HBD.
    HBD: nitrogen-hydrogen (N-A) and oxygen-hydrogen (O-H) groups capable of donating a hydrogen bond.


    Rule: compound should have no more than 10 HBA.
    HBA: N and O atoms cabable of accepting an H-bond.

    Consequences of more H-bonds: 
        - Increased Solubility (Good): H will bond with water molecules.
        - Poor pPermeability (Bad): too many H-bonds prevent ability to cross cell membranes, reducing bioavailability.
        - Decreased Oral Bioavailability (Bad): not easily absorbed in gastrointestinal tract
    Vice versa for less H-bonds.

    Tuned to prefer 2 HBD and 3.6 HBA, and penalize breaking of rule. 
    """
    def __init__(self, Wandb = True):
        super(LipinskiReward,self).__init__(Wandb)
        
    def name(self, Wandb = True):
        return 'LipinskiReward'

    def rescale(self, x):
        # keep in range [-2, 2]
        return x
    
    def _reward(self, mol):
        numHdonors = Lipinski.NumHDonors(mol)
        numHacceptors = Lipinski.NumHAcceptors(mol)
        # print(f"(numHdonors, numHacceptors) = {(numHdonors, numHacceptors)}")

        penalty = 0
        # penalize for breaking lipinskis rule of 5
        if numHdonors > 5:
            penalty += 1
        if numHacceptors > 10:
            penalty +=1
        # reward for being similar to "best" molecules
        targetHdonors = 2
        HBD_dist = 1
        targetHacceptors = 6
        HBA_dist = 4

        # targetHdonors = 2
        # HBD_dist = 2.5
        # targetHacceptors = 5
        # HBA_dist = 4

        reward = SingleReward.gaussianReward(numHdonors, targetHdonors, HBD_dist) + SingleReward.gaussianReward(numHacceptors, targetHacceptors, HBA_dist)
        # return - penalty
        return reward - penalty

class QEDReward(SingleReward):
    """
    Quantitative Evaluation of Drug-Likeness (QED) reward evaluates drug-likeness based on 
        - Lipophilicity (LogP)
        - Molecular Weight (MW)
        - HBD and HBA
        - Polar Surface Area (PSA): Surface area of molecule that is polar
        - # Rotatable Bonds: Flexibility of molecule. 
        - # Aromatic Rings: Presence of aromaticity (eg. Benzene: 6 carbon atoms w/ alternating double and single bonds)

    Raw reward range: [0, 1)

    Tuned to prefer > 0.5
    """
    def __init__(self, Wandb=True):
        super(QEDReward,self).__init__(Wandb)
        
    def name(self):
        return 'QED Reward'
    
    def rescale(self, x):
        # emphasize vals above 0.5
        transf_r = SingleReward.transform(x, steepness=8)
        # translate to range [-2,2] from [0,1]
        return SingleReward.scale_to_range(transf_r)
        # return SingleReward.scale_to_range(transf_r, new_range=(-4,4))
    
    def _reward(self, mol):
        try: 
            r = qed(mol)
            # print(f"qed(mol) = {qed(mol)}")
            return r
        except:
            return 0

class MolBuffer():
    """Priority Queue for keeping track of best molecules seen so far
    """

    def __init__(self, length: int):
        """_summary_

        Args:
            length (int): number of molecules to hold
        """
        self.length = length
        self.__buffer = []
        self.step = 1

        # wandb.define_metric("best_mols",step_metric="best_mols_step")

    def put(self, item):
        """Put (score,smile) into buffer. Will not replace anything if its bad

        Args:
            item (tuple): score, smile tuple
        """
        if len(self.__buffer) < self.length:
            heapq.heappush(self.__buffer, item)
        else:
            heapq.heappushpop(self.__buffer, item)

    def log(self, top: int):
        """log top molecules to wandb run

        Args:
            top (int): how many molecules to log
        """
        top_n = heapq.nlargest(top,self.__buffer, lambda x: x[0])
        for item in top_n:
            try:
                molecule = Chem.MolFromSmiles(item[1])
                pil_image = Draw.MolToImage(molecule, size=(300, 300))            
                wandb.log({'best_mols': wandb.Image(pil_image)})
            except:
                pass
            

    def reset(self):
        self.__buffer = []  
    
class FinalRewardModule():
    '''reward function handling'''
    def __init__(self,writer,r_list,wandb_log = True, scaling=False):
        self.r_list = r_list
        self.writer = writer
        self.buffer = MolBuffer(1)
        self.wandb_log = wandb_log

        self.color = ['blue','orange','red','green','yellow', 'pink', 'purple']
        
        self.n_iter = 0
        self.scaling = scaling
        self.gamma = .99
        self.mean_discount = .95
        self.mean = 0
        self.first = True
        """Maybe put a buffer of rewards and use that to calculate std
        """
        
        if scaling:
           self.r_stats = [0]
    
    def UpdateTrainingModules(self):
        pass
    
    def print_rewards(self):
        for idx, reward in enumerate(self.r_list):
            print(f"Reward {idx + 1} - {reward.name()}: {reward.reward_list}")
            
    def PlotRewards(self):
        for idx,SingleReward in enumerate(self.r_list):
            plt.plot(SingleReward.reward_list, label = SingleReward.name(),color = self.color[idx])
        plt.show()
        plt.savefig('rewards_plot.png')  # Save the plot as an image
        
    
 


    def LogRewards(self):
        self.buffer.log(5)
        self.buffer.reset()
          
    def GiveReward(self,mol):
        if self.wandb_log:
            wandb.log({'num_mols':self.n_iter})
        self.n_iter += 1
        mol.UpdatePropertyCache()
        rewards = 0
        for rewardObject in self.r_list:
            # print(f"giving reward type {rewardObject} to mol {Chem.MolToSmiles(mol)}")
            reward = rewardObject.giveReward(mol)
            # self.writer.add_scalar(rewardObject.name(), reward, self.n_iter)
            rewards += reward

        self.buffer.put((rewards,Chem.MolToSmiles(mol)))
        
        if self.scaling:         
            if self.first:
                self.mean = 0
                self.first = False
                
            else:
                self.mean = rewards*(1-self.mean_discount) + (self.mean* self.mean_discount)
            
            rewards_centered = (rewards-self.mean)#/(np.std(self.r_stats) + 1e-10)
            if self.wandb_log:
                wandb.log({'running norm rewards': rewards_centered})
            return rewards_centered
        if self.wandb_log:
            wandb.log({'total rewards': rewards})
        return rewards
            
            
            
def main():
    mol = Chem.MolFromSmiles('Cc1c2c3c(c(NPCl)c1PC2)CN=N3')
    sr = SynthReward(False)
    print(sr.giveReward(mol))

    fr = FinalRewardModule()

if __name__ == '__main__':
    main()








