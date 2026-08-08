#!/bin/bash
#SBATCH --job-name=gudi_MaCh3Input
#SBATCH --output=slurm_files/GUDI_log/gudi_MaCh3Input-%A_%a.out
#SBATCH --error=slurm_files/GUDI_log/gudi_MaCh3Input-%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Flynn.Y.Guo@stonybrook.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --array=1-48    # one task per file
#SBATCH --time=02:00:00

# Load any required modules (if necessary)
# On SeaWulf
#  module load root

 # On NNhome machine
source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh

# Adjust for the local data location.
export OA_INPUTS=/storage/shared/DUNE/OA-inputs
export OA_OUTPUTS=/storage/shared/fyguo/Work/Projects/GUNDAM/gudi_v4/

# Declare an (initially empty) array
files=()

# Loop through the splits manually
for file in ${OA_INPUTS}/MaCH3-inputs/v2/DUNE_2023_FD_CAFs/FD_*; do
    files+=("$file")
done


# Use the SLURM_ARRAY_TASK_ID to select the file from the array.
selected_file=${files[$((SLURM_ARRAY_TASK_ID - 1))]}

# Run your program and pass the selected file as an argument.
echo "./macros/prepareGundamMCTree_MaCh3Input $selected_file "
./macros/prepareGundamMCTree_MaCh3Input "$selected_file"
