#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00 
#SBATCH --mem=50G

# at seawulf
source /gpfs/home/uyevarouskay/Work/env_gundam_v211_v1.sh #gundam version 2.1.1 (installed 2026-04-17)

# at SBU NN home cluster
#source /home/isgould/work/gundam/install/setup.sh

# at seawulf cluster
export OA_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v5/

# at dunegpvm cluster
#export OA_INPUTS=/pnfs/dune/persistent/users/uyevarou/GUNDAM/OA-inputs/lbl/

# at SBU NN home cluster
#export OA_INPUTS=/storage/shared/DUNE/OA-inputs/lbl/gudi-inputs/v5/

export NUOSCILLATOR_ROOT_LIB=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh


# if you request an Asimov fit, simply disregard the `-d` option, --scan is to build llh scans
gundamFitter -a -c ./config_DUNE.yaml -t 8
# To run few files only consider override the parameters config:
#gundamFitter -a -d -c ./config_DUNE.yaml -t 8 -of overrides/fewfiles.yaml
# To disable  oscillation parameters consider override the parameters config:
#gundamFitter -a -d --scan -c ./config_DUNE.yaml -of overrides/disableOscillationParameters.yaml -t 8

