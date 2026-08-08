#!/bin/bash
# Update the "external" shared libraries that are used by dials to
# calculate weights and load events.  This can be safely run multiple
# times.  This must be run after gundam has been configured.
#
# Note: This is where the NuOscillator interface is compiled.
#

if ! which gundamFitter > /dev/null; then
    echo GUNDAM not properly configured and gundamFitter not found.
    echo Externals for gudi-anoa-sens25a not updated.
    exit 1
fi

if ! which root-config > /dev/null; then
    echo ROOT not properly configured and root-config not found.
    echo Externals for gudi-anoa-sens25a not updated.
    exit 1
fi

# Adjust if necessary.  This could be done more intelligently, but
# should be OK.
CXX=$(root-config --cxx)
CXX_FLAGS=$(root-config --cflags)
CXX_FLAGS+=" -shared"
CXX_FLAGS+=" -fPIC"

# Not required for DUNE
# ${CXX} ${CXX_FLAGS} ./inputs/datasets/variables/pmuCC.cpp -o ./inputs/datasets/variables/pmuCC.so
# ${CXX} ${CXX_FLAGS} ./inputs/datasets/variables/invalidateDuplicateEntry.cpp -o ./inputs/datasets/variables/invalidateDuplicateEntry.so

# Synchronize and recompile the oscillator interface.
git submodule update
(cd gundamOscAnaTools/resources/TabulateNuOscillator; bash compile.sh)
