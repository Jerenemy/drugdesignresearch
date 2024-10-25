# Check if the package is installed and, if not, install it to a custom directory
if (!require(ggplot2)) install.packages("ggplot2", lib="~/R/x86_64-pc-linux-gnu-library")
if (!require(tidyr)) install.packages("tidyr", lib="~/R/x86_64-pc-linux-gnu-library")
if (!require(gridExtra)) install.packages("gridExtra", lib="~/R/x86_64-pc-linux-gnu-library")

# Load the packages from the custom directory
library(ggplot2, lib.loc="~/R/x86_64-pc-linux-gnu-library")
library(tidyr, lib.loc="~/R/x86_64-pc-linux-gnu-library")
library(gridExtra, lib.loc="~/R/x86_64-pc-linux-gnu-library")

# Read the data from the CSV file

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a CSV file path was provided
if (length(args) == 0) {
  stop("Error: No CSV file path provided. Please provide the CSV file as an argument.")
}

# The first argument is the CSV file path
csv_file <- args[1]

# Check if the CSV file exists
if (!file.exists(csv_file)) {
  stop(paste("Error: CSV file", csv_file, "not found!"))
}


data <- read.csv(csv_file)
# data <- read.csv('/home/jeremy/RFiles/csvs/bad_rewards3.csv')
# data <- read.csv("/home/jeremy/RFiles/csvs/control_real_mols.csv")
# data <- read.csv('/home/jeremy/RFiles/csvs/chembl_25.csv')

# data <- read.csv("/home/jeremy/RFiles/csvs/bad_rewards.csv")
# data <- read.csv("/home/jeremy/RFiles/csvs/bad_rewards1.csv")
# data <- read.csv("/home/jeremy/RFiles/csvs/real_mol_rewards.csv")
# data <- read.csv('/home/jeremy/RFiles/csvs/only_pk11000.csv')


#data <- read.csv("/home/jeremy/RFiles/csvs/molecules_rewards.csv")


# Set the Nickname column as a factor with levels in the order they appear in the CSV file
data$Nickname <- factor(data$Nickname, levels = data$Nickname)

# Dynamically get the list of reward types to plot separately
reward_types <- grep("Reward", names(data), value = TRUE)

# Define a set of colors
colors <- c("blue", "red", "green", "purple", "orange", "darkgreen")

# Specify the directory to save the plots
save_dir <- "../Files"

# Create the directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Create a list to store the plots
plots <- list()

# Loop through each reward type and create a separate plot
for (i in seq_along(reward_types)) {
  reward <- reward_types[i]
  color <- colors[i %% length(colors) + 1]
  
  # Create the plot
  p <- 
    ggplot(data, aes_string(x = "Nickname", y = reward)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 0, alpha = 0.2, fill = "red") +  # UPDATE: Add red translucent background below zero
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2, alpha = 0.2, fill = "green") + # UPDATE: Add green translucent background above zero
  
  
    geom_bar(stat = "identity", fill = color) +
    labs(title = paste("Reward:", reward), x = "Molecule", y = "Reward Amount") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),  # Reduce text size
          plot.margin = margin(1, 1, 2, 1, "cm")) +  # Add margin to avoid text overlap
    ylim(-2, 2)  # UPDATE: Set y-axis limits to ±2 and center at 0
    
  # Add the plot to the list
  plots[[i]] <- p
}

# Arrange all plots in a grid
combined_plot <- do.call(grid.arrange, c(plots, ncol = 2))

# Save the combined plot as a PNG file

ggsave(filename = paste(save_dir, '/PRESENTATION_affinity_rewards.png', sep = ''), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/bad_rewards3.png', sep = ''), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/control_real_mols.png', sep = ''), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/chembl_25.png', sep = ''), plot = combined_plot, width = 20, height = 20)

# ggsave(filename = paste(save_dir, '/bad_rewards.png', sep = ''), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/bad_rewards1.png', sep = ''), plot = combined_plot, width = 20, height = 20)


#ggsave(filename = paste(save_dir, "/all_in_1_plot_good.png", sep = ""), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, "/all_in_1_plot_bad.png", sep = ""), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, "/all_in_1_plot_bad1.png", sep = ""), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/all_in_1_plot_complx.png', sep = ""), plot = combined_plot, width = 20, height = 20)
# ggsave(filename = paste(save_dir, '/only_pk11000.png', sep = ''), plot = combined_plot, width = 20, height = 20)

