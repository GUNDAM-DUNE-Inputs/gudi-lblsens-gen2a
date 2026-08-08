#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00 
#SBATCH --mem=50G

# at seawulf
source /gpfs/home/uyevarouskay/Work/env_gundam_v211.sh #gundam version 2.1.1 (installed 2026-02-04)

# inputs
export OA_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v4_joint/
export ND_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/Work/ND_workshop/

# at seawulf cluster
export NUOSCILLATOR_ROOT_LIB=/gpfs/projects/McGrewGroup/uyevarou/Work/atm/gudi-anoa-sen25a/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

# if you request an Asimov fit, simply disregard the `-d` option, --scan is to build llh scans
gundamFitter -a -d -c ./config_DUNE_ND.yaml -t 8 #--light-mode

