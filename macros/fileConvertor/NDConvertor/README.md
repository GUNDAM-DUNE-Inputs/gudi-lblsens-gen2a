# Running the ND Flat-CAF-to-GUNDAM Converter at FNAL

This package provides the source code for the
`GUNDAMNDLArAnalysis_flatCAF` converter.

The converter reads DUNE ND flat CAF files and produces GUNDAM-compatible
inputs for the ND-LAr analysis.

## Contents

This directory contains:

```text
README.md
GUNDAMNDLArAnalysis_flatCAF.cxx
```

The converter source must be copied into a local checkout of
[`DUNE/duneanasel`](https://github.com/DUNE/duneanasel) before building.

> **Important**
>
> `GUNDAMNDLArAnalysis_flatCAF.cxx` is not currently included in the public
> `DUNE/duneanasel` repository. In addition to copying the source file, the
> user must add the corresponding CMake target to
> `examples/ndlarplot/CMakeLists.txt`.

## Prerequisites

You need:

- A Fermilab account with DUNE access.
- Access to the required input files in `/pnfs/dune`.
- Apptainer available through CVMFS.
- A local checkout of `DUNE/duneanasel`.
- A local checkout of `DUNE/ND_CAFMaker`.

The commands below use:

```bash
WORKDIR=/exp/dune/app/users/<username>
```

Replace `<username>` with your Fermilab username.

> **Note**
>
> The input-file paths below point to a specific MiniProdN5p3 FHC production.
> Access to these files may depend on your DUNE and PNFS permissions.

## 1. Clone the required repositories

```bash
cd /exp/dune/app/users/<username>

git clone https://github.com/DUNE/duneanasel.git
git clone https://github.com/DUNE/ND_CAFMaker.git
```

`ND_CAFMaker` provides `ndcaf_setup_deps.sh`, which configures the UPS-based
runtime dependencies used to build and run the converter.

## 2. Start an SL7 Apptainer environment

Run the following from your FNAL working directory:

```bash
cd /exp/dune/app/users/<username>

apptainer exec \
  --bind /cvmfs:/cvmfs \
  --bind /pnfs:/pnfs \
  --bind /exp:/exp \
  --bind $PWD:$PWD \
  --pwd $PWD \
  /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest \
  bash
```

The remaining commands must be run inside this container.

## 3. Set up ND CAFMaker dependencies

```bash
source /exp/dune/app/users/<username>/ND_CAFMaker/ndcaf_setup_deps.sh prof
```

Use `prof` for an optimized build. Use `debug` instead when debugging:

```bash
source /exp/dune/app/users/<username>/ND_CAFMaker/ndcaf_setup_deps.sh debug
```

## 4. Copy the converter source

Copy `GUNDAMNDLArAnalysis_flatCAF.cxx` from this package into the
`ndlarplot` example directory:

```bash
cp /path/to/GUNDAMNDLArAnalysis_flatCAF.cxx \
  /exp/dune/app/users/<username>/duneanasel/examples/ndlarplot/
```

Confirm that the file was copied successfully:

```bash
ls -lh /exp/dune/app/users/<username>/duneanasel/examples/ndlarplot/GUNDAMNDLArAnalysis_flatCAF.cxx
```

## 5. Register the converter in CMake

The converter must be registered as a CMake executable target.

Open:

```text
/exp/dune/app/users/<username>/duneanasel/examples/ndlarplot/CMakeLists.txt
```

and add the CMake block used in the development version of the converter.

> **Important**
>
> The exact CMake block depends on the libraries linked by the converter.
> Copy the `GUNDAMNDLArAnalysis_flatCAF` target definition from the development
> checkout where the converter was originally built. A source file by itself
> cannot be compiled unless CMake is instructed to create and link the target.

After the change, verify that the target name appears in the file:

```bash
grep -n "GUNDAMNDLArAnalysis_flatCAF" \
  /exp/dune/app/users/<username>/duneanasel/examples/ndlarplot/CMakeLists.txt
```
## 6. Configure and build the converter

From the `duneanasel` checkout, configure the build directory and build the
flat-CAF converter:

```bash
cd /exp/dune/app/users/<username>/duneanasel

cmake -S . -B build \
  -Dduneanasel_USE_UPS_duneanaobj=ON \
  -Dduneanasel_build_examples=ON

# Build only the flat-CAF-to-GUNDAM converter.
cmake --build build --target GUNDAMNDLArAnalysis_flatCAF -j 4
```

The executable should be located at:

```text
/exp/dune/app/users/<username>/duneanasel/build/examples/ndlarplot/GUNDAMNDLArAnalysis_flatCAF
```

Check that it was built:

```bash
ls -lh build/examples/ndlarplot/GUNDAMNDLArAnalysis_flatCAF
```

## 7. Convert one flat CAF file

```bash
CONVERTER=/exp/dune/app/users/<username>/duneanasel/build/examples/ndlarplot/GUNDAMNDLArAnalysis_flatCAF

INPUT_FILE=/pnfs/dune/persistent/users/jlarkin/ndprod/MiniProdN5/\
MiniProdN5p3_NDComplex_FHC.nusyst.full.sanddrift.spineonly.v2/\
CAF.flat/0000000/\
MiniProdN5p3_NDComplex_FHC.nusyst.full.sanddrift.spineonly.v2.0000248.cafnusyst.flat.root

"$CONVERTER" "$INPUT_FILE"
```

Check the converter output and configuration to determine the location and name
of the produced GUNDAM input file.

## 8. Convert all flat CAF files in a directory

Run the converter separately for each input file:

```bash
CONVERTER=/exp/dune/app/users/<username>/duneanasel/build/examples/ndlarplot/GUNDAMNDLArAnalysis_flatCAF

INPUT_DIR=/pnfs/dune/persistent/users/jlarkin/ndprod/MiniProdN5/\
MiniProdN5p3_NDComplex_FHC.nusyst.full.sanddrift.spineonly.v2/\
CAF.flat/0000000

for input_file in "$INPUT_DIR"/*.cafnusyst.flat.root; do
  [ -e "$input_file" ] || continue
  "$CONVERTER" "$input_file"
done
```

Do not pass a shell wildcard directly to the converter unless it has been
verified that the executable supports multiple input-file arguments. The loop
above processes one file at a time and is safer for large productions.

## Notes

- This workflow is intended for interactive testing and modest-scale
  conversion.
- For large productions, submit one job per CAF file or a small group of CAF
  files through the FNAL batch system.
- Record the `duneanasel` commit, `ND_CAFMaker` commit, converter version, and
  input production used for each GUNDAM production.
- A Git branch containing the source and its CMake configuration is preferred
  for long-term reproducibility.