# Running macros

## 1. Prepare Inputs for GUNDAM
- ```prepareGundamMCTree_MaCh3Input.cpp```: read files from MaCh3's input

This repository contains a script that prepares files in a GUNDAM-friendly format by copying the `cafTree` and adding new branches, including the converted dials into a `TClonesArray` of `TGraphs`.

### New Branches Added
The following branches are added to the `event_tree`:

- `POTScaledWeight`
- `POT_generated`
<!-- - `POTweight` -->
- `Nonswap`
- `Nueswap`
- `Tauswap`
- `isNC`
- `isRHC`
- `sample_idx`
- `Converted dials`

### Compiling the Script

To compile the `prepareGundamMCTree_MaCh3Input.cpp` script, use the following command:

```bash
g++ -o prepareGundamMCTree_MaCh3Input prepareGundamMCTree_MaCh3Input.cpp -I$(root-config --incdir) $(root-config --libs) -std=c++17
```

This command compiles the C++ code and links it against ROOT libraries. Make sure ROOT is properly installed and configured on your system.

### Running the Script Locally

Once the script is compiled, use the `run_jobs.sh` script to execute it. Adjust the data location in the script.

```bash
chmod +x run_jobs.sh
./run_jobs.sh
```

### Submitting Jobs to SLURM: running the Slurm Array Job

1. **Edit the Slurm script**  

   - Update the `#SBATCH` directives (partition, time, memory, CPUs) to match your cluster’s available resources.  
   - Set the `OA_INPUTS` and `OA_OUTPUTS` variables to point at your data directories.

2. **Submit the job**  

```bash
sbatch run_jobs_slurm.sh
```
