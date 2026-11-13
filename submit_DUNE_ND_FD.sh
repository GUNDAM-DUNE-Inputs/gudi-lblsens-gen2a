#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00 
#SBATCH --mem=50G

# at seawulf
source /gpfs/projects/McGrewGroup/uyevarou/env_gundam_v211_v2.sh #gundam version v2.1.1f+7-g8324995c/main (installed 2026-05-14)

# inputs
export OA_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v4_joint/
#export OA_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v5/
export ND_INPUTS=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/ND_MiniProd53_v1/

# at seawulf cluster
export NUOSCILLATOR_ROOT=/gpfs/projects/McGrewGroup/uyevarou/Work/atm/gudi-anoa-sen25a/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT}/bin/setup.NuOscillator.sh

# if you request an Asimov fit, simply disregard the `-d` option, --scan is to build llh scans
gundamFitter -a -d -c ./config_DUNE_ND_FD.yaml -t 8 #--light-mode

