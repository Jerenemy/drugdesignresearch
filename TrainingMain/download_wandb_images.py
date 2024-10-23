import wandb
import os
raise Exception("don't run this way too many images to download")

# Authenticate with wandb
api = wandb.Api()

# Define your project and run details
entity = "jzay-Wesleyan University"  # Your wandb username or team name
project = "DrugNet"  # Your project name
run_id = "3eysqqba"  # Your run ID

# Create a directory to save images
os.makedirs("wandb_images", exist_ok=True)

# Get the run
run = api.run(f"{entity}/{project}/{run_id}")

# Get all the files from the run
files = run.files()

# Loop through the files and download images
for file in files:
    if file.name.endswith(".png") or file.name.endswith(".jpg"):
        file.download(root="wandb_images", replace=True)

print("All images downloaded!")
