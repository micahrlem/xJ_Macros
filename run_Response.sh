#!/bin/bash

# Script to run make_xJ_Response_smear with various JES/JER configs and input samples

MACRO="make_xJ_Response_smear_flatten"
RVAL=0.2

# JES/JER config files
CONFIGS=( "binning_posJES.config" "binning_negJES.config" )

# Input ROOT files
INFILES=(
   "pythia8-Jet10-Run21-multiR.root"
   "pythia8-Jet20-Run21-multiR.root"
   "pythia8-Jet30-Run21-multiR.root"
   # "Herwig-Jet10-Run21-multiR.root"
    # "Herwig-Jet30-Run21-multiR.root"
)

# Loop through all combinations
for CONFIG in "${CONFIGS[@]}"; do
    for INFILE in "${INFILES[@]}"; do
        echo "Running macro with config: $CONFIG and infile: $INFILE"
        root -l -b -q "${MACRO}.C(\"${CONFIG}\", ${RVAL}, \"${INFILE}\")"
    done
done
