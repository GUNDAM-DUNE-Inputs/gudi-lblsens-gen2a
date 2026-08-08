# GUNDAM Inputs and Configuration for DUNE Long-baseline Neutrino OA
## 0. Input files
* **On NN home:**

  ```bash
  export OA_INPUT_FOLDER=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v5/
  ```

* **On dunegpvm:**

  ```bash
  export OA_INPUT_FOLDER=/pnfs/dune/persistent/users/uyevarou/GUNDAM/OA-inputs/lbl/gudi-inputs/v5/
  ```

* **On seawulf:**

  ```bash
  export OA_INPUT_FOLDER=/gpfs/projects/McGrewGroup/uyevarou/gudi-inputs/lbl/v5/
  ```

## 1. Set Up the GUNDAM Environment ⚙️


To install GUNDAM, follow [GUNDAM Getting Started](https://gundam-organization.github.io/gundam/GettingStarted.html).

Then, in `submit_DUNE.sh`, source your installation. Replace with your actual install path:

```bash
source yourpath/gundam/install/setup.sh
```

For users on the NN home cluster you can source pre-installed GUNDAM version 2.X:

```bash
source /home/fyguo/GUNDAM/Install/gundam/setup.sh
```

Alternatively, you can add gundam to your path as described in [GUNDAM Getting Started](https://gundam-organization.github.io/gundam/GettingStarted.html):

```bash
export PATH="$INSTALL_DIR/gundam/bin:$PATH"
export LD_LIBRARY_PATH="$INSTALL_DIR/gundam/lib:$LD_LIBRARY_PATH"
```

## 2. Compile gundamOscAnaTools🔧

Run the following command:

```bash
# Clone the main repo
git clone https://github.com/GUNDAM-DUNE-Inputs/gudi-lbl-tdr25a
cd gudi-lbl-tdr25a

# Initialize and update submodules
git submodule init
git submodule update
# Compile
./update-externals.sh
```
Alternatively you can build manually. If needed, configure `./gundamOscAnaTools/resources/TabulateNuOscillator/CMakeLists.txt` and compile by using:

```bash
# before compiling source gundam
source yourpath/gundam/install/setup.sh
# if there are no GPUs modify ./gundamOscAnaTools/resources/TabulateNuOscillator/CMakeLists.txt replace:
#-set(UseProbGPULinear 1 CACHE STRING  "NuOscillator Option -- Enable ProbGPULinear")
#+set(UseProbGPULinear 0 CACHE STRING  "NuOscillator Option -- Enable ProbGPULinear")
# if the compilation fails add: 
#+target_link_libraries(dumpAtmOscProb PRIVATE ${CMAKE_DL_LIBS})
#+target_link_libraries(dumpLblOscProb PRIVATE ${CMAKE_DL_LIBS})
#+target_link_libraries(testTabulatedNuOscillator PRIVATE ${CMAKE_DL_LIBS})

cd ./gundamOscAnaTools/resources/TabulateNuOscillator/
mkdir build-$(uname -m)
cd build-$(uname -m)
cmake -DCMAKE_INSTALL_PREFIX=${PWD} ..
make install
```


## 3. Use NuOscillator Interface (instructions for submodule) 🔧

To add `gundamOscAnaTools` to your path make sure that the following lines appear in `submit_DUNE.sh`:

```bash
export NUOSCILLATOR_ROOT=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT}/bin/setup.NuOscillator.sh
```

## 4. Running GUNDAM

Use correctly configured `submit_DUNE.sh` script to perform either a dry run (```-d```) or a full Asimov fit (```-a```). You can adjust the number of cores according to your cluster using the `-t` option.

Submit the job script:
* **On NN home:**
```bash
sbatch nnhome_slurm.sh
```

## 5. Enable/disable/modify Oscillation Parameters ⚙️

**Disable Oscillation Parameters** 

  The Oscillation parameters are provided in the `./inputs/parameters/configParSet.yaml`. To run GUNDAM without oscillation parameters you can use an override file by running the following command

  ```bash
  gundamConfigUnfolder -c config_DUNE.yaml -of overrides/disableOscillationParameters.yaml -o
  ```

  It will create an override configuration file `config_DUNE_With_disableOscillationParameters.json` that can be run by modifying `submit_DUNE.sh` in the following way:

  ```bash
  gundamFitter -a -d -c ./config_DUNE_With_disableOscillationParameters.json -t 8
  ```

  Alternatively you can disable oscillation parameter set in the `./inputs/parameters/configParSet.yaml` by changing the following flag to `false` (not recommended):

  ```yaml
  - name: "Oscillation Parameters"
    isEnabled: false
  ``` 
  And run `submit_DUNE.sh` without modification.

  ```bash
  ./submit_DUNE.sh
  ```

**Change the priors** 


  Different priors can be provided for the oscillation parameters by using corresponding override files. For instance, to set up prior values to PDG recommended use the following override file:

  ```bash
  gundamConfigUnfolder -c config_DUNE.yaml -of overrides/pdg25NormalOscillationParameters.yaml -o
  ```
 
  Then modify `submit_DUNE.sh` to:

  ```bash
  gundamFitter -a -d -c ./config_DUNE_With_pdg25NormalOscillationParameters.json -t 8
  ```

  And run:

  ```bash
  ./submit_DUNE.sh
  ``` 

---

## Using NuOscillator Interface (old instructions for not a submodule case) 🔧

To use the GUNDAM interface for NuOscillator, follow these steps:

**Build the NuOscillator Interface**
Make sure that GUNDAM is set up and in your PATH. You will also need to have the development environment set up. Then navigate to the directory and compile the interface:

```bash
cd gundamOscAnaTools/resources/TabulateNuOscillator/
./compile.sh
```

**On NN home:**

```bash
source /home/uyevarou/Work/env_dune_lbl.sh
cd gundamOscAnaTools/resources/TabulateNuOscillator/
```

**On dunegpvm:**

```bash
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12
spack load cmake@3.27.7%gcc@12.2.0
spack load gcc@12.2.0
# First time compile only to properly get ROOT exports
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/bin/thisroot.sh
```

**Modify the GUNDAM source file and add additional libraries to your PATH** (included in `env_dune_lbl.sh`):

```bash
source /YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/bin/setup-NuOscillator.sh
# or alternatively:
export LD_LIBRARY_PATH="/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/_deps/oscprob-src/lib/:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/lib/:$LD_LIBRARY_PATH"
```

---

**Modify the Run Line** 📝

In the `submit_DUNE.sh` script, update the run line as follows:

```bash
export NUOSCILLATOR_ROOT=/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/
export NUOSCILLATOR_ROOT_LIB=/YOURPATH/gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
```

---

## Run GUNDAM 🎉

```bash
./submit_DUNE.sh
```

Happy fitting! 😊
