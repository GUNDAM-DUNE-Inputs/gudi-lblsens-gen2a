#!/usr/bin/env bash
#SBATCH -p long-40core-shared
#SBATCH --array=0-41                           # 6 params × (2·NUM_SIGMA+1)= 42;  (NUM_SIGMA = 2) σ-points ⇒ 30 tasks
#SBATCH --cpus-per-task=8
#SBATCH -t 4:00:00
#SBATCH --mem=50G
set -euo pipefail
set -x

# ── Config ────────────────────────────────────────────────────────────────────
input_yaml="./overrides/Profiling.yaml"
output_dir="/gpfs/projects/McGrewGroup/uyevarou/outputs/DUNE/lbl/profiling/output_frozen_TDR/"

# parameter definitions (must stay in sync)
paramName=(sin12 sin13 sin23 m21 m32 dcp)
occurrence=(1    2     3     4   5    6  )
priorDefault=(0.310 0.0224 0.582 7.39E-5 0.002525 0.0)
sigma=(        0.013  0.0007   0.0217 1.8e-6     3.4e-5     1.046 ) #- 7 points

# NUM_SIGMA=10 only dcp
#paramName=(dcp)
#occurrence=( 6  )
#priorDefault=( 0 )
#sigma=(  0.314  ) #1.047

NUM_SIGMA=3                                #  ±3σ ⇒ 7 points (-2σ,-1σ,0,+1σ,+2σ)
#NUM_SIGMA=10                              # ±10σ ⇒ 21 points
STEP_COUNT=$(( 2*NUM_SIGMA + 1 ))               # =5
PARAM_COUNT=${#paramName[@]}                   # =6

# ── Map SLURM task → (param index, σ-offset) ─────────────────────────────────
TASK_ID=$SLURM_ARRAY_TASK_ID                    # 0…29
param_idx=$(( TASK_ID / STEP_COUNT ))           # integer division
step_idx=$(( TASK_ID % STEP_COUNT ))            # 0…4
offset=$(( step_idx - NUM_SIGMA ))              # maps 0→−2,1→−1,2→0,3→+1,4→+2

# ── Pull out the right values ────────────────────────────────────────────────
name=${paramName[param_idx]}
occ=${occurrence[param_idx]}
default=${priorDefault[param_idx]}
step=${sigma[param_idx]}

# ── Compute the new prior = default + offset·σ ───────────────────────────────
newPrior=$(awk -v d=$default -v s=$step -v k=$offset \
  'BEGIN{ printf "%g\n", d + k*s }')

echo "occ=$occ, prior=$newPrior, name=$name"

echo "Task#$TASK_ID → $name (occ=$occ), offset=$offset σ ⇒ prior=$newPrior"
#echo "./submit_DUNE_profiling.sh "$occ" "$newPrior" "$name" "
# ── Launch your profiling script ─────────────────────────────────────────────
./submit_DUNE_profiling.sh "$occ" "$newPrior" "$name" "$input_yaml" "$output_dir"