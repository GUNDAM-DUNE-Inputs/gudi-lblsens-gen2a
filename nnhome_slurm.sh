#!/bin/bash
#SBATCH --job-name=GundamInputDUNE_LBL
#SBATCH --output=slurm_files/GundamInput_log/GundamInputDUNE_LBL_%j.out
#SBATCH --error=slurm_files/GundamInput_log/GundamInputDUNE_LBL_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=Flynn.Y.Guo@stonybrook.edu
#SBATCH --mem=50G
#SBATCH -c 16
#SBATCH --gres=gpu

source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh
source /home/fyguo/GUNDAM/Install/gundam/setup.sh

#display detailed information about your system
uname -a
#display the date
date
echo "started"

export OA_INPUTS=/storage/shared/DUNE/OA-inputs/lbl/gudi-inputs/v5

export LD_LIBRARY_PATH="/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/_deps/oscprob-src/lib/:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/lib/:$LD_LIBRARY_PATH"

export NUOSCILLATOR_ROOT=/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64
export NUOSCILLATOR_ROOT_LIB=/home/fyguo/Work/Projects/GUNDAM/GundamInputDUNE_LBL/Devel/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64

# run fitter
# 1. run profiling:
# echo "gundamFitter -d -a -c config_DUNE_With_PDG_oscparameter_flat.json --inject-toy-parameter ./injector/parameterInjector1.yaml -t 16 --toy 0"
# gundamFitter -d -a -c config_DUNE_With_PDG_oscparameter_flat.json --inject-toy-parameter ./injector/parameterInjector1.yaml -t 16 --toy 0
# 2. run asimov
echo "gundamFitter -c config_DUNE.yaml -a -d --gpu -t 16"
gundamFitter -c config_DUNE.yaml -a -d --gpu -t 16

echo "finished"
date