# Generating CSVs of Molecules with Rewards

This file contains instructions for Generating CSVs of Molecules with Rewards in the format I have been using.

## External Links

- [Wandb](https://wandb.ai/jzay/projects)
- [Important Notes](https://docs.google.com/document/d/1J4tggblSIhgjd74qx43tRQQQVuRp6-jsfC5kksxItus/edit?usp=sharing)
- [Plots](https://drive.google.com/drive/folders/193gJ95C34S1ehT_ABo4qHrAFoppiMLZ-?usp=sharing)
- [ChatGPT](https://chatgpt.com/)
- [Smiles to Image](http://hulab.rxnfinder.org/smi2img/)

## Table of Contents

- [Generating CSVs of Molecules with Rewards](#generating-csvs-of-molecules-with-rewards)
  - [External Links](#external-links)
  - [Table of Contents](#table-of-contents)
  - [Setup](#setup)
  - [CSV Creation](#csv-creation)
  - [Rewards](#rewards)
  - [Images Creation](#images-creation)
  - [Graphing with R](#graphing-with-r)
  - [Overall Improvements](#overall-improvements)

## Setup
**SSH**:
ssh into Theo's pc (if not using it).
- ssh jzay@129.133.102.210
- Lob$ter66


For now, can only ssh from another PC with a wired connection. What if I connect the PC to eduroam? Then can I ssh in from my mac
- I connected to temp_devices (password: wesdevice2023 (i think))
- Now can ssh into PC from my mac when connected to temp_devices, but not sure if it works with the vpn (if im off campus)

**Conda Environment:** Must use Theo's conda env *prawn*.
- /home/theo/anaconda3/condabin/conda init
- conda activate prawn

## CSV Creation

**CWD:** When creating formatted CSVs, must do it from the *Jeremy/CSV* folder:
- cd /home/jzay/Desktop/GoodThesis1/DrugDesignThesis/CLEAN/Jeremy/CSV

**[create_csv.sh](Jeremy/CSV/create_csv.sh)**: Create your CSV file by running a bash script -  for now this does NOT take cmd-line args, edit the arg line(s) in the script.
- ./create_csv.sh

create_csv.sh first cd's into the *GenAITutorial/Training/TrainingData* directory, where it calls either [get_lines.py](../../../../GenAITutorial/Training/TrainingData/get_lines.py) or [get_random_lines.py](../../../../GenAITutorial/Training/TrainingData/get_random_lines.py).

**[get_lines.py](../../../GenAITutorial/Training/TrainingData/get_lines.py)**: selects specific lines from an input CSV file and writes them to an output CSV file. Does not change the formatting.
This file needs to be edited to be in the Jeremy/CSV folder.

**[get_random_lines.py](../../../GenAITutorial/Training/TrainingData/get_random_lines.py)**: same thing as get_lines.py but random lines instead of specific lines.

create_csv.sh then copies the created CSV file back to *Jeremy/CSV/files* and cd's back into *Jeremy/CSV*.

create_csv.sh then calls [convert_csv_to_list.py](CSV/convert_csv_to_list.py).

**[convert_csv_to_list.py](Jeremy/CSV/convert_csv_to_list.py)**: the bash script passes the output_csv and the file_type as cmd-line args, and this script writes relevant data from the CSV into 2 vars as a var name and a tuple. 
- What is the purpose of this step?
  - To have a list that can be rerun to calculate the rewards, without having to recreate the csv file?
  - Is this step unneccessary?
    - Could be skipped: calculate the reward straight from the new csv file created. 
    - Bad because every time a new list is created it adds a new one to the file, even if it is a duplicate, adding unneccessary clutter.
  - Instead, have func that returns this variable, call this func when creating new csv file. 
  - Output new csv file instead of overwriting.
  - Have separate folders: one for input CSVs, one for output CSVs. 

**[csv_driver.py](Jeremy/CSV/csv_driver.py)**: bash script calls csv_driver.py and passes the output_csv created in the previous step as a cmd-line arg. csv_driver.py updates the rewards fields in the csv file by calling funcs in [rewards.py](../Rewards/rewards.py) for each line (molecule). 

## Rewards

**[rewards.py](Rewards/rewards.py)**: goal is to tweak rewards to make each reward more discriminatory and useable. 
- Each reward still needs to be scaled correctly by having the most important rewards weighted more (increase their range in the rescale)
- Have a few functions to work with to construct the rewards:
  - **[gaussianReward](Rewards/rewards.py#L58)** calculates reward of *x* from gaussian dist, based on *ideal* and *width*.
  - **[scale_to_range](Rewards/rewards.py#L80)** scales *x* from orig finite range to new finite range.
  - **[transform](Rewards/rewards.py#L105)** sigmoid tranformation enphasizes vals above/below *threshold* based on *steepness*.
    - Expected range is [0,1].
  - **[normalize_value](Rewards/rewards.py#L135)** scales *x* from orig range to [0,1]. 
    - Although this is the exact same method as scale_to_range, it is necessary to be called before passing into *transform*.
  - **[combinedReward](Rewards/rewards.py#L157)** calculates reward starting as linear up to *transition_point*, then smoothly transitions to gaussian centered around *ideal* val with fwhm *sigma*. ]
    - The point of combinedReward is to give enough incentive to allow the model to build up to close to that point.
    - With pure gaussian reward the vals would be too small at first and never climb, as the differences would be insignificant.


**[SynthReward](Rewards/rewards.py#L178)**: calculates synthetic ability (how easy is it to create) of the molecule.
- Raw reward: calculates in range [1,9] (1 is best).
- Scaled reward: normalized, transformed to emphasize above 0.37, scaled to range [-2,2]


**[DockReward](Rewards/rewards.py#L226)**: calculates docking score between molecule and y220c.
- Raw reward: Decent is around -5 (lower is better).
  - If calculation takes more than 5 seconds, end process and give penalty (+3).
- Scaled reward: can't scale since no finite range, instead take negative, subtract 5, clip to range [-2,10].


**[SizeReward](Rewards/rewards.py#L280)**: calculates reward based on size of molecule.
- Raw reward: returns num atoms in molecule.
- Scaled reward: use *combinedReward*, since doesn't stay as really small. 
  - target_atoms = 24
  - dist = 10
  - transition_point = 20
  - scale_to range [-2,2] from [0,1]


**[LogPReward](Rewards/rewards.py#L315)**: calculates how hydrophilic/lipophilic the molecule is.
- Raw reward: centered around 0
  - LogP < 0 ==> more hydrophobic (lipophilic): less soluble, breaks through cell membranes better.
  - more hydrophilic (lipophobic): more soluble, doesn't break through cell membranes as well.
- Scaled reward: *gaussianReward* target 3.3 (tuned to prefer hydrophobic). Scaled to range [-2,2].
  

**[LipinskiReward](Rewards/rewards.py#L346)**: 
- Compound should have no more than 5 HBD.
  - HBD: nitrogen-hydrogen (N-A) and oxygen-hydrogen (O-H) groups capable of donating a hydrogen bond. 
- Compound should have no more than 10 HBA.
  - HBA: N and O atoms cabable of accepting an H-bond.
- Consequences of more H-bonds: 
    - Increased Solubility (Good): H will bond with water molecules.
    - Poor pPermeability (Bad): too many H-bonds prevent ability to cross cell membranes, reducing bioavailability.
    - Decreased Oral Bioavailability (Bad): not easily absorbed in gastrointestinal tract
- Vice versa for less H-bonds.
- Raw reward: reward - penalty
  - reward: gaussianReward centered around target HBD/HBA (sum of both) (+1 for being correct for each)
  - penalty: breaking the rule (-1 for each infraction)
- Scaled reward: already in range [-2,-2], because maximum reward = +2, minimum penalty = -2.


## Images Creation

**[imgs_driver.py](Jeremy/CSV/imgs_driver.py)**: [create_csv.sh](CSV/create_csv.sh) calls imgs_driver.py with the output_csv as a cmd-line arg, and displays the first 25 mols in the CSV on one page, 5 mols per line. Outputs the img [here](CSV/Output/ai_output/control_real_mols1.png) with name f"{output_csv}{i}".
- Could be improved by taking cmd-line args to specify tot num mols graphed (num mols per line, num lines) and resolution.


## Graphing with R

Instructions: 
1. Download CSV file to R CSV folder.
2. Locate R script to run to produce desired graph.
3. Copy 2 lines of printed output to terminal from running csv_driver.py and replace lines in R script to read CSV file and save graph.
4. Drag graph produced to google drive. 

Could be improved by having bash script take output of csv_driver.py (lines to edit in R script), pass to R script as cmd-line args, edit R script, run R script to produce graph. 
- Would require running R script on Theo's computer (while ssh'd in).
- Would be easier if could run R script in VSCode so can see clearly while ssh'd in.
  - Not neccessary, especially when running directly off of Theo's computer.


## Overall Improvements

1. Fix cd'ing in bash script, shouldn't be neccessary.
2. Fix storing list as var in separate file.
3. Weigh rewards correctly.
4. Improve/replace AutoDock Vina bad docking calculation. 
5. Make bash script to automate graphing with R.
6. Obtain output of list of best smiles strings created by model (not just imgs from WandB).
7. Make fingerprint matrix reward that evaluates reward based on *features*, how many features molecule has similar to target molecule. 
8. Run MD on best molecule to see if it's actually docking.

How to add location-specific docking information back to model? Problem I'm having is that the reward model just takes a value as the reward. I can't decide what specifically it should change. 

Figure out what PPO is, what step-reward is, anything else I don't understand yet.

Went through [sv_utils.py](TrainingUtils/SupervisedTraining/sv_utils.py), kinda understand it. still need to understand how changing the args affects it. 

Still need to understand exactly  the order of the stages. What is called when i run it with [main.py](TrainingMain/main.py)? main.py calls [config.yaml](TrainingMain/config_small.yaml). Where is [make_session.py](TrainingMain/make_session.py) called? What about [training.py](TrainingMain/training.py)? [PPO.py](TrainingUtils/PPO.py)?

## How Theo's Code Works
Idk yet, I'll find out tomorrow.

## How to Fix to make Reward Creation Much More Streamlined
Should be able to do everything by putting a CSV file into the INPUT folder, running the command with the csv_name as the cmd-line arg, the specific task to perform (eg img creation, reward evaluation, R graph creation), and it should spit the CSV with rewards in the OUTPUT, along with the Images in the OUTPUT, and any graphs created in the OUTPUT. This shouldn't be too difficult.

## Updates made for Thayer presentation
I attempted to update the generation process, but I was up against the clock and at the last minute my generation process wasn't actually running the docking reward, was just giving -2 each time. So I went back to the old version and was able to calculate the rewards for the new effectors. This was like 2 weeks ago or more though, and I can barely remember.

