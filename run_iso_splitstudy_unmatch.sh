#!/bin/bash

# Define generators and jet types
generators=("pythia8" "Herwig")
jets=("Jet10" "Jet20" "Jet30")
Rs=("0.2" "0.4")

# Loop over all combinations
for gen in "${generators[@]}"; do
  for jet in "${jets[@]}"; do

    # Skip invalid HERWIG Jet20 case
    if [[ "$gen" == "Herwig" && "$jet" == "Jet20" ]]; then
      continue
    fi

    for R in "${Rs[@]}"; do
      infile="${gen}-${jet}-Run21-multiR.root"
      echo "Running for: ${infile}, R = ${R}"
      root -l -b -q "Make_Iso_Split_Study_UnmatchTruth.C(\"${infile}\", \"binning_original.config\", ${R})"
    done

  done
done
